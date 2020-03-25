// graphics.cc - Simple windows
// $Id: graph.cc,v 1.1 2012/03/08 19:32:10 jack Exp jack $

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;

#include "graph.h"
#include "icon.cc"

#define	EVENT_BITS	(KeyPressMask | StructureNotifyMask | \
			 ButtonPressMask | ExposureMask)

#define	GC_BITS		(GCFunction   | GCForeground | \
   			 GCFillStyle  | GCArcMode)

#define	WIN_BITS	(CWBackPixel  | CWBorderPixel)

#define	BORDER_MAX		10
#define	BORDER_MIN		0

#define	A1	  		0
#define	A2			23040


#define	DEF_BACKGROUND		"White"
#define	DEF_BORDERCOLOR		"Black"
#define	DEF_FOREGROUND  	"Black"

#define	DEF_OUTLINE		1
#define	DEF_TRACE		0

#define	DEF_BORDER		5
#define	DEF_GEOMETRY		"=400x400-50+50"

Display*  Xcon::display;
int	  Xcon::screen;
int	  Xcon::depth;
Window	  Xcon::rootWindow;
Colormap  Xcon::colormap; 

int Xcon::count = 0;
int Xcon::active = 0;
char *Xcon::app_name = NULL;

Xcon::Xcon(char* name, int activate)
{
  count++;
  app_name = name;
  if (activate && !active) connect();
}

Xcon::Xcon(int activate)
{
  count++;
  if (activate && !active) connect();
}

int Xcon::connect() {
   /****************
    * Get defaults *
    ****************/
  if (active) return 1;
  display = XOpenDisplay(NULL);
  if (display == NULL) {
    cerr << app_name << ": Cannot open display" << endl;
    active = 0;
    exit( 1 );
  }
  active = 1;

  screen     = XDefaultScreen(display);
  rootWindow = XDefaultRootWindow(display);
  depth      = XDefaultDepth(display, screen);
  if (depth > 1) colormap = XDefaultColormap( display, screen );

  return active;
}

void Xcon::disconnect()
{
  if (active) XCloseDisplay(display);
  active = 0;
}

char* Xcon::Default(const char* par) const
{
   return XGetDefault(display, app_name, par);
}

/********************
 * Allocate a color *
 ********************/

long Xcon::AllocColor(char* colorName, char* cDefault) const
{
  XColor	xcolor, dummy;
  Status	status;
  long color;
  if (colorName == 0) colorName = "";
  if (depth > 1) {
    status = XAllocNamedColor(display, colormap, colorName,
			      &xcolor, &dummy);
    if (status)
      color = xcolor.pixel;
    else {
      status = XAllocNamedColor(display, colormap, cDefault,
				&xcolor, &dummy);
      if (status)
	color = xcolor.pixel;
      else
	color = -1;
    }
    if (color < 0)
      color = BlackPixel(display, screen);
  }
  else
    color = BlackPixel(display, screen);
  return color;
}

long Xcon::AllocCol(int col) const
{
  XColor	xcolor, dummy;
  Status	status;
  long color;
  xcolor.red = col*10; xcolor.green = col*10; xcolor.blue = 20;
  status = XAllocColor(display, colormap, &xcolor);
  if (status)
    color = xcolor.pixel;
  else
    color = BlackPixel(display, screen);
  if (color == 0 && col > 0) color = WhitePixel(display, screen);
  return color;
}

/******************************
 * Initialize the color table *
 ******************************/

void Xcon::InitColorTable(char* backgr, char* foregr, char* border) {
  if (depth > 1) {
    background = AllocColor(backgr, DEF_BACKGROUND);
    if (background < 0)
      background = WhitePixel( display, screen );
    foreground = AllocColor(foregr, DEF_FOREGROUND);
    if (foreground < 0)
      foreground = BlackPixel( display, screen );
    bordercoul = AllocColor(border, DEF_BORDERCOLOR);
    if (bordercoul < 0)
      bordercoul = BlackPixel( display, screen );
  }
  else {
    background  = WhitePixel( display, screen );
    foreground  = BlackPixel( display, screen );
    bordercoul  = BlackPixel( display, screen );
  }
}
Xcon::~Xcon() { if (--count == 0) disconnect(); }

Win::Win(char* name, int activate) : Xcon(name,activate) { open(); }

Win::Win() { open(); }

Win::~Win() { close(); }

int Win::open()
{
    if (!active) return 0;

    /*******************
     * Allocate colors *
     *******************/
    
    char * backgr = Default("background");
    char * foregr = Default("foreground");
    char * border = Default("borderColor");
    InitColorTable(backgr, foregr, border);

    /*****************************
     * Create the window *
     *****************************/

    char * geometry     = Default("geometry" );

      /**********************************
       * Allow window to only be square *
       **********************************/

    XSizeHints hints;
    //    hints.min_aspect.x = hints.min_aspect.y = 1;
    //    hints.max_aspect = hints.min_aspect;
    hints.flags=0; //USPosition; // PAspect

    bwidth = 0;
    XWMGeometry(display, screen, geometry, DEF_GEOMETRY,
		bwidth,  &hints, &x, &y, &w, &h, &hints.win_gravity);

    XSetWindowAttributes	setWinAtt;
    setWinAtt.background_pixel = background;
    setWinAtt.border_pixel     = bordercoul;

    window = XCreateWindow(display, rootWindow, x, y, w, h,
			   bwidth,  depth, InputOutput,
			   CopyFromParent, WIN_BITS, &setWinAtt);

    Pixmap icon_pixmap =	// Define icon.
      XCreateBitmapFromData(display,window,icon_bits,icon_width,icon_height);
    XWMHints wm_hints;
    wm_hints.initial_state=NormalState;
    wm_hints.input=True;
    wm_hints.icon_pixmap=icon_pixmap;
    wm_hints.flags=StateHint|IconPixmapHint|InputHint;

    XSetWMProperties(display,window,0,0,0,0, &hints, &wm_hints,0);
    XStoreName(display, window, app_name);

    /* Map Window and wait for it to apear... */
    XSelectInput(display, window, ExposureMask);
    XMapWindow(display, window);
    XWindowEvent(display, window, ExposureMask, &event);
    XSync(display, 0);

    /*********************
     * Adjust the window *
     *********************/

    XWindowAttributes		winAtt;
    XGCValues			gcValues;

    XGetWindowAttributes(display, window, &winAtt);

    w = winAtt.width;
    h = winAtt.height;

    /****************************************
     * Create and set the Graphics contexts *
     ****************************************/

    gcValues.function   = GXcopy;
    gcValues.foreground = foreground;
    gcValues.fill_style = FillSolid;
    gcValues.arc_mode   = ArcPieSlice;

    gc = XCreateGC(display, window, GC_BITS, &gcValues);

    drawable = window;		// Graphics to window by default.
				// Alternative is to a pixmap as follows:
    // Create pixmap //
    // drawable = pixmapW = XCreatePixmap(display, window, w, h, depth);
    // XFillRectangle(display, pixmapW, gcN, 0, 0, w, h);

    // Select input //
    event_mask = EVENT_BITS;
    XSelectInput(display, window, EVENT_BITS);

    XFlush(display);

    // Flags //
    mapped = 1;

    xscale = yscale = 1;
    xmargin = ymargin = 0;

    return 1;
}

void Win::close()
{
//  XFreePixmap(display, pixmapW);
  if (!active) return;
  XDestroyWindow(display, window);
  XFreeGC(display,gc);
  mapped = 0;
  XFlush(display);
}

void Win::Clear() {
  XClearWindow(display, window);
}

void Win::Redraw() {}

int Win::SetFont(char *name, GC gct) {
  XFontStruct *font = XLoadQueryFont(display,name);
  if (!font) font = XLoadQueryFont(display,"-*-helvetica-bold-r-normal--8-80-75-75-p-50-iso8859-1");
  if (!font) font = XLoadQueryFont(display,"8x16");
  if (!font) {
    cerr << app_name <<
      " ::Win : Cannot open font " << name << endl;
    exit(-1);
  }
  if (gct) XSetFont(display,gct,font->fid);
  else     XSetFont(display,gc ,font->fid);
  return font->ascent+font->descent;
}

int Win::DrawString(int x, int y, const char *s, GC gc) {
  if (gc == 0) gc = Win::gc;
  return XDrawString(display,drawable,gc, x, y, s,strlen(s));
}

// ProcessEvents looks for user events (i.e. keyboard presses) so
// that it can toggle display flags (or quit).

void Win::ProcessEvents()
{
  KeySym	keysym;
  char		buff[1];

  while (XPending(display) &&
	 XCheckWindowEvent(display, window, event_mask, &event)) {

    switch (event.type ) {
    case (KeyPress):
      XLookupString(&event.xkey, buff, 1, &keysym, NULL );
      switch (buff[0]) {
      case 'o':
      case 'O':
	break;
      case 'Q':
	exit(0);
      case 'l':
      case 'L':
	break;
      case 'a':
      case 'A':
	break;
      case 'c':
      case 'C':
	Clear();
	break;
      case 'w':
      case 'W':
	break;
      case 'S':
      case 's':
	break;
      case 'F':
	break;
      case 'f':
// 	  DrawCircle(event.xkey.x-2,event.xkey.y-2,4);
// 	  XCopyArea(display, pixmapW, window, gcN, 0, 0, w, h, 0, 0);
// 	  XSync(display, 0);
// 	  curves.deltaF(x,y);
	//cout << " deltaF = "
	//     << curves.closest(x,y)->deltaF(0.001,yr/200) << endl;
	break;
      default:
	break;
      }
      break;

    case (ConfigureNotify):

      // Window is being resized.

      w    = event.xconfigure.width;
      h    = event.xconfigure.height;
      break;

    case (ButtonPress):
      break;

    case (Expose):
      if (event.xexpose.count == 0) {
	Redraw();
      }
      break;

    default:
      break;
    }
  }
}

void Win::DrawCircle(int x, int y, int size, GC gc) const
{
  // Draw Circle.
  if (gc==0) gc = Win::gc;
  x -= size/2; y -= size/2;
  XDrawArc(display, drawable, gc, x, y, size, size, A1, A2 );
  if (x + size > w)
    XDrawArc(display, drawable, gc, x-w, y, size, size, A1, A2 );
  if (y + size > h)
    XDrawArc(display, drawable, gc, x, y-h, size, size, A1, A2 );
}

void Win::DrawSolidCircle(int x, int y, int size, GC gc) const
{
  // Draw the filled in circle.
   if (gc==0) gc = Win::gc;
   x -= size/2; y -= size/2;
   XFillArc(display, drawable, gc, x, y, size, size, A1, A2 );
   if (x + size > w)
     XFillArc(display, drawable, gc, x-w, y, size, size, A1, A2 );
   if (y + size > h)
     XFillArc(display, drawable, gc, x, y-h, size, size, A1, A2 );
}

void Win::DrawRectangle(int x, int y, int w, int h, GC gc) const
{
  if (gc==0) gc = Win::gc;
  // Draw the filled in rectangle.
  XFillRectangle(display, drawable, gc, x, y, w, h);
}

void Win::plot(double x, double y, GC gc)
{
  if (gc==0) gc = Win::gc;
  XDrawPoint(display, drawable, gc, xc(x), yc(y));
}

void Win::line(double x1, double y1, double x2, double y2, GC gc)
{
  if (gc==0) gc = Win::gc;
  XDrawLine(display, drawable, gc, xc(x1), yc(y1), xc(x2), yc(y2));
}

void Win::arrow(double x1, double y1, double x2, double y2, GC gc)
{
  line(x1,y1,x2,y2,gc);
  double dx = x2-x1, dy = y2-y1, size = 1; // sqrt(dx*dx+dy*dy);
  double vx = 0.4*dx - dy*0.2, vy = 0.4*dy + dx*0.2;
  line(x2,y2, x2-size*vx,y2-size*vy,gc);
  vx = 0.4*dx + dy*0.2; vy = 0.4*dy - dx*0.2;
  line(x2,y2, x2-size*vx,y2-size*vy,gc);
}

void Win::warrow(double x1, double y1, double x2, double y2, GC gc)
{
  Drawable temp = drawable;
  drawable = window;
  line(x1,y1,x2,y2,gc);
  double dx = x2-x1, dy = y2-y1, size = 1; // sqrt(dx*dx+dy*dy);
  double vx = 0.4*dx - dy*0.2, vy = 0.4*dy + dx*0.2;
  line(x2,y2, x2-size*vx,y2-size*vy,gc);
  vx = 0.4*dx + dy*0.2; vy = 0.4*dy - dx*0.2;
  line(x2,y2, x2-size*vx,y2-size*vy,gc);
  drawable = temp;
}

void Win::plot(double *x, double *y, int n,
	       double x0, double y0, double x1, double y1,
	       GC gc)
{
  if (gc==0) gc = Win::gc;
  double lx = x1 - x0, ly = y1 - y0;
  XPoint *points = new XPoint[n];
  for (int i = 0; i < n; i++) {
    points[i].x = short(w*(x[i] - x0)/lx);
    points[i].y = short(h - h*(y[i] - y0)/ly);
  }
  XDrawLines(display, drawable, gc, points, n, CoordModeOrigin);
  delete [] points;
}

void Win::circle(double x, double y, double size, int solid, GC gc)
{
  if (solid)
    DrawSolidCircle(xc(x), yc(y), int(size*w*xscale), gc);
  else
    DrawCircle(xc(x), yc(y), int(size*w*xscale), gc);
}

void Win::waitformouse()
{
  XRaiseWindow(display,window);
  key = keysym = mouse = 0;
  do {
    XWindowEvent(display,window, ButtonPressMask, &event);
  } while (!( event.type == ButtonPress &&
	      event.xbutton.window == window ));
  mouse = event.xbutton.button;
  mouse_x = event.xbutton.x;
  mouse_y = event.xbutton.y;
}

Window Win::SubWin(int x,int y,unsigned int w,unsigned int h,
		   unsigned int bwidth,
		   unsigned long bordercoul, unsigned long background) {
  return XCreateSimpleWindow(display, window, x, y, w, h,
			     bwidth, bordercoul, background);
}

TextWin::TextWin() : cx(0),cy(16),fh(10){
  XSetWindowAttributes	setWinAtt;
  setWinAtt.backing_store = WhenMapped;
  XChangeWindowAttributes(display, window, CWBackingStore, &setWinAtt);
}

void TextWin::Text(char *text, GC gc)
{
  DrawString(cx, cy, text, gc);
  cy += fh;
}

void VortexWin::Info(int rader, char *s) {
  if (info==0) return;
  int cx = fh, cy = 2*fh;
  for (; rader > 0; rader--) {
    char *sl = strchr(s,'\n');
    XDrawString(display,drawable,gc, cx, cy, s,sl-s);
    cy += fh; s = sl+1;
  }
}

VortexWin::VortexWin() : howf(20), showc(0), please_close(0), no3DL(0) {
  XGCValues gcValues;
  gcValues.function   = GXcopy;
  gcValues.foreground = background;
  gcValues.fill_style = FillSolid;
  gcValues.arc_mode   = ArcPieSlice;
  gcB = XCreateGC(display, window, GC_BITS, &gcValues);

  gcValues.foreground = AllocColor("red","black");
  gcbox = XCreateGC(display, window, GC_BITS, &gcValues);
  XSetLineAttributes(display, gcbox, 6, LineSolid, CapRound, JoinRound);

  gcValues.foreground = AllocColor("red","black");
  gcred = XCreateGC(display, window, GC_BITS, &gcValues);
  XSetLineAttributes(display, gcred, 2, LineSolid, CapRound, JoinRound);

  gcValues.foreground = AllocColor("black","black");
  gc3 = XCreateGC(display, window, GC_BITS, &gcValues);
  XSetLineAttributes(display, gc3, 4, LineSolid, CapRound, JoinRound);
 
  gcValues.foreground = AllocColor("blue","black");
  gcdis = XCreateGC(display, window, GC_BITS, &gcValues);
  XSetLineAttributes(display, gcdis, 4, LineSolid, CapRound, JoinRound);
 
  gcblue = XCreateGC(display, window, GC_BITS, &gcValues);

  gcValues.foreground = AllocColor("green","black");
  gcgreen = XCreateGC(display, window, GC_BITS, &gcValues);

  // Create pixmap //
  drawable = XCreatePixmap(display, window, w, h, depth);
  XFillRectangle(display, drawable, gcB, 0, 0, w, h);

  info = 0; topview = 0; translate_back = 0; show_connections = 0;
  fh = SetFont("-adobe-new century schoolbook-bold-i-normal--18-180-75-75-p-111-iso8859-1");

  xmm = ymm = 1;
  xm = ym = 2;
  repl = -1;
  replmax = 0;

  xmargin = ymargin = 0.5;
}

VortexWin::~VortexWin() { if (active) XFreePixmap(display, drawable); }

void VortexWin::ShowTemp(double T, char *str) {
  if (!active) return;
  char s[40];
  if (str) sprintf(s,"%s: T = %lg (%s)",app_name,T,str);
  else     sprintf(s,"%s: T = %lg",app_name,T);
  XStoreName(display, window, s);
}

void VortexWin::ShowTitle(const char *str) {
  if (!active) return;
  char s[60];
  if (str) sprintf(s,"%s: ( %s )",app_name,str);
  else     sprintf(s,"%s",app_name);
  XStoreName(display, window, s);
}

int VortexWin::Events()
{
  if (XPending(display) &&
      XCheckWindowEvent(display, window,
  			KeyPressMask | ButtonPressMask, &event)) {

    key = keysym = mouse = 0;

    switch (event.type ) {
    case (KeyPress):
      XLookupString(&event.xkey, &key, 1, &keysym, NULL );
      switch (key) {
      case 'l':
      case 'L':
	no3DL = 1-no3DL;
	key = 0;
	break;
      case '>':
	xp += 0.5; break;
      case '<':
	xp -= 0.5; break;
      case 'R':
	zp += 0.5; break;
      case 'r':
	zp -= 0.5; break;
      case 'Z':
	yp += 0.5; break;
      case 'z':
	yp -= 0.5; break;
      case 'Q':
	exit(0);
      case 'X':
	please_close = 1; break;
      case 'i':
      case 'I':
	info = 1 - info;
	break;
      case 't':
      case 'T':
	topview = (1+topview) % 3;
	break;
      case 'M':
	xmm *= 1.2; ymm *= 1.2; break;
      case 'm':
	xmm /= 1.2; ymm /= 1.2; break;
      case 'u':
	howf = 2*howf; break;
      case 'n':
	howf = (howf > 1) ? howf/2 : 1; break;
      case 'w':
	show_connections = (1 + show_connections) % 4;
	break;
      case 'B':
      case 'b':
	translate_back = 1 - translate_back;
	break;
      default:
	break;
      }
      break;

    case (ButtonPress):
      mouse = event.xbutton.button;
      mouse_x = event.xbutton.x;
      mouse_y = event.xbutton.y;
      break;

    default:
      break;
    }
    return 1;
  }
  return 0;
}

void VortexWin::Draw()
{
  while (XPending(display) &&
	 XCheckWindowEvent(display, window,
			   StructureNotifyMask | ExposureMask, &event)) {

    switch (event.type ) {

    case (ConfigureNotify):

      // Window is being resized.

      w    = event.xconfigure.width;
      h    = event.xconfigure.height;
      XFreePixmap(display, drawable);
      drawable = XCreatePixmap(display, window, w, h, depth );
      XFillRectangle(display, drawable, gcB, 0, 0, w, h);
      break;

    case (Expose):
      if (event.xexpose.count == 0) {
	Win::Redraw();
 	XCopyArea(display, drawable, window, gc, 0, 0, w, h, 0, 0);
 	XSync(display, 0);
      }
      break;

    default:
      break;
    }
  }
  // Show...
  XCopyArea(display, drawable, window, gc, 0, 0, w, h, 0, 0);
  // Clear...
  XFillRectangle(display, drawable, gcB, 0, 0, w, h);
}

void VortexWin::ProcessEvents()
{
  KeySym	keysym;
  char		buff[1];

  while (XPending(display) &&
	 XCheckWindowEvent(display, window, event_mask, &event)) {

    switch (event.type ) {
    case (KeyPress):
      XLookupString(&event.xkey, buff, 1, &keysym, NULL );
      switch (buff[0]) {
      case 'u':
      case 'U':
	howf = 2*howf;
	break;
      case 'n':
      case 'N':
	howf = (howf > 1) ? howf/2 : 1;
	break;
      case 'o':
      case 'O':
	no3DL = 1-no3DL;
	break;
      case '>':
	xp += 0.5; break;
      case '<':
	xp -= 0.5; break;
      case 'R':
	zp += 0.5; break;
      case 'r':
	zp -= 0.5; break;
      case 'Z':
	yp += 0.5; break;
      case 'z':
	yp -= 0.5; break;
      case 'Q':
	exit(0);
      case 'X':
	please_close = 1;
	break;
      case 'i':
      case 'I':
	info = 1 - info;
	break;
      case 'c':
      case 'C':
	Clear();
	break;
      case 't':
      case 'T':
	topview = 1 - topview;
	break;
      case 'w':
      case 'W':
	show_connections = 1 - show_connections;
	break;
      case 'M':
	xmm *= 1.2; ymm *= 1.2; break;
      case 'm':
	xmm /= 1.2; ymm /= 1.2; break;
      case 'B':
      case 'b':
	translate_back = 1 - translate_back;
	break;
      case '1':
	repl++;
	if (repl > replmax) repl = 0;
	break;
      case '2':
	repl--;
	if (repl < 0) repl = replmax;
	break;
      case '3':
	repl = -1;
	break;
      case 'S':
      case 's':
	break;
      case 'F':
	break;
      case 'f':
// 	  DrawCircle(event.xkey.x-2,event.xkey.y-2,4);
// 	  XCopyArea(display, pixmapW, window, gcN, 0, 0, w, h, 0, 0);
// 	  XSync(display, 0);
// 	  curves.deltaF(x,y);
	//cout << " deltaF = "
	//     << curves.closest(x,y)->deltaF(0.001,yr/200) << endl;
	break;
      default:
	break;
      }
      break;

    case (ConfigureNotify):

      // Window is being resized.

      w    = event.xconfigure.width;
      h    = event.xconfigure.height;
      XFreePixmap(display, drawable);
      drawable = XCreatePixmap(display, window, w, h, depth );
      XFillRectangle(display, drawable, gcB, 0, 0, w, h);
      break;

    case (ButtonPress):
//       if (event.xbutton.button == 1) {
// 	curves.ymax *= 1 - double(event.xbutton.y)/h;
// 	curves.xmax *= double(event.xbutton.x)/w;
// 	Show();
//       }
//       else if (event.xbutton.button == 2) {
// 	curves.ymax /= 1 - double(event.xbutton.y)/h;
// 	curves.xmax /= double(event.xbutton.x)/w;
// 	Show();
//       }	
      // mouse at (event.x, event.y)
      // quit = 1;
      break;

    case (Expose):
      if (event.xexpose.count == 0) {
	Redraw();
 	XCopyArea(display, drawable, window, gc, 0, 0, w, h, 0, 0);
 	XSync(display, 0);
      }
      break;

    default:
      break;
    }
  }
  // Show...
  XCopyArea(display, drawable, window, gc, 0, 0, w, h, 0, 0);
  // Clear...
  XFillRectangle(display, drawable, gcB, 0, 0, w, h);
}

 void VortexWin::setview(double xmi, double ymi, double xma, double yma) {
  xmin = xmi;
  ymin = ymi;
  xmax = xma;
  ymax = yma;
  xscale = 1/(xmax-xmin);
  yscale = 1/(ymax-ymin);
}

void VortexWin::setview(double L) {
  xp = L/3; yp = +1; zp = +L/3; lll = L;
  xscale = yscale = 1/L;
}

void VortexWin::tr3D(double x, double y, double z, double &xw, double &yw) {
  if (topview) {
    double yy = yp; if (yy <= 0) yy = 0.1;
    xw = x/yy; yw = y/yy; return;
  }
  y*=0.05;
  x -= xp; y += yp; z-= zp;
  if (y <= 0) { xw=0; yw=0; return; }
  x /= y; z/= y; x+= xp; y-=yp; z+=zp;
  xw = x; yw = z;
}
void VortexWin::tr3Dp(double x, double y, double z, int &xi, int &yi) {
  double xx,yy; tr3D(x,y,z,xx,yy);
  xi = xc(xx/xm); yi = yc(yy/ym);
}

void VortexWin::dr3Dbox(double L, double Z) {
  double xx[16] = { 0, L, L, 0, 0, 0, L, L, L, L, L, L, 0, 0, 0, 0 };
  double yy[16] = { 0, 0, L, L, 0, 0, 0, 0, 0, L, L, L, L, L, L, 0 };
  double zz[16] = { 0, 0, 0, 0, 0, Z, Z, 0, Z, Z, 0, Z, Z, 0, Z, Z };
  double xxx[16],yyy[16]; xm=ym=0;
  for (int i = 0; i < 16; i++) {
    tr3D(xx[i],yy[i],zz[i],xxx[i],yyy[i]);
    if (xxx[i] > xm) xm = xxx[i];
    if (yyy[i] > ym) ym = yyy[i];
  }
  xm *= xmm; ym *= ymm; //??
  XSetLineAttributes(display, gcred, 2,
		     LineSolid, CapRound, JoinRound);
  plot(xxx,yyy,16,0,0,xm,ym,gcred);
}

void VortexWin::plot3D(double x1,double y1,double z1,
		       double x2,double y2,double z2, GC gc3) {
  if (gc3==0) gc3 = Win::gc;
  int xx1, yy1, xx2, yy2;
  tr3Dp(x1,y1,z1,xx1,yy1);
  tr3Dp(x2,y2,z2,xx2,yy2);
  if (!no3DL) {
    int line_width = int( 10/((y1+y2)*0.1+yp)+1 );
    if (topview) line_width = xc(0.25/lll/yp);
    if (line_width < 0) line_width = 0;
    XSetLineAttributes(display, gc3, line_width,
		       LineSolid, CapRound, JoinRound);
  }
  XDrawLine(display,drawable,gc3,xx1,yy1,xx2,yy2);
}

void VortexWin::plot3Dsz(double x1,double y1,double z1,
			 double x2,double y2,double z2,
			 double size, GC gc3) {
  if (gc3==0) gc3 = Win::gc;
  int xx1, yy1, xx2, yy2;
  tr3Dp(x1,y1,z1,xx1,yy1);
  tr3Dp(x2,y2,z2,xx2,yy2);
  if (!no3DL) {
    int line_width = int( size*10/((y1+y2)*0.1+yp)+1 );
    if (topview) line_width = int( size*xc(0.25/lll/yp) );
    if (line_width < 0) line_width = 0;
    XSetLineAttributes(display, gc3, line_width,
		       LineSolid, CapRound, JoinRound);
  }
  XDrawLine(display,drawable,gc3,xx1,yy1,xx2,yy2);
}
