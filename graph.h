// graph.h - Simple windows
// $Id: graph.h,v 1.1 2003/03/25 15:09:02 jack Exp $

#ifndef Graph_H
#define Graph_H

#include <X11/Xlib.h>
#include <X11/Xutil.h>

// A Win is the object that manages all X11 display activities.

// Handles the X-connection

class Xcon {
public:
  // Static variables
  static char*	  app_name;
  static Display  *display;
  static int	  screen;
  static int	  depth;
  static Window	  rootWindow;
  static Colormap colormap; 

  long background, foreground, bordercoul;

  static int	  active;
  // private:
  static int      count;	// How many instances?

public:
  Xcon() { if (count++ == 0 || !active) connect(); }
  virtual ~Xcon();

  Xcon(char *name, int activate = 1);
  Xcon(int activate);

  int connect();
  void disconnect();

  char*		Default(const char*) const;
  long		AllocColor(char*, char*) const;
  long		AllocCol(int col) const;
  void		InitColorTable(char* backgr, char* forgr, char* border);
}; 

class Win : public Xcon {
public:

  Window	window;
  Drawable      drawable;	// Window or pixmap.
  XEvent	event;
  long          event_mask;

  int		mapped;

  int		x, y, w, h;
  int		bwidth;

  int		mouse, mouse_x, mouse_y;
  char		key;
  KeySym	keysym;

  GC		gc;

  double	xscale, yscale;	// Scale factor
  double	xmargin, ymargin; // shift window

  Win(char *name, int activate = 1);
  Win();
  virtual ~Win();

  int	 open();
  void	 close();

  void   Clear();
  virtual void   Redraw();

  int    SetFont(char *name, GC gc = 0); // Returns hight of chars.
  int    DrawString(int x, int y, const char *s, GC gc = 0);

  void   ProcessEvents();
  void	 DrawCircle(int, int, int, GC gc = 0) const;
  void	 DrawSolidCircle(int, int, int, GC gc = 0) const;
  void	 DrawRectangle(int, int, int, int, GC gc = 0) const;

  // Translate coords to window coords.
  int    xc(double xx) { return (int) ((xx+xmargin)*w*xscale); }
  int    yc(double yy) { return h - (int) ((yy+ymargin)*h*yscale); }

  // Translate back.
  double xc(int xx) { return double(xx)/(w*xscale)-xmargin; }
  double yc(int yy) { return double(h - yy)/(w*yscale)-ymargin; }

  void   plot(double  x, double  y, GC gc = 0); // Plot one point.
  void   line(double  x, double  y, double x2, double y2, GC gc = 0);
  void   arrow(double  x, double  y, double x2, double y2, GC gc = 0);
  void   warrow(double x, double  y, double x2, double y2, GC gc = 0);
  void   plot(double *x, double *y, int n,
	      double x0, double y0, double x1, double y1,
	      GC gc = 0); // Plot several points. xi is corners.
  void   circle(double x, double y, double size, int solid = 0, GC gc = 0);
  void   waitformouse();
  void   Flush() { XFlush(display); }
  Window SubWin(int x,int y,unsigned int w,unsigned int h,
		unsigned int bwidth,
		unsigned long bordercoul, unsigned long background);
};

class TextWin : public Win {
public:
  int cx, cy, fh;			// current position, font height.
  TextWin();
  virtual ~TextWin() {}
  int SetFont(char *name, GC gcc = 0) { return fh = Win:: SetFont(name,gcc); }
  void Text(char *text, GC gc = 0);
//  void Redraw();
};

class VortexWin : public Win {
public:
  VortexWin();
  virtual ~VortexWin();

  GC  gcB, gcbox, gcdis;
  GC  gcred, gcblue, gcgreen, gc3;

  int howf,showc;
  int please_close, no3DL, translate_back, repl, replmax;

  int timetoshow(int id)
  { return (active && (repl == -1 || id == repl) &&
	    (++showc >= howf ? showc = 0, 1 : 0)); }

  int timetoshow() { return (active && (++showc >= howf ? showc = 0, 1 : 0)); }

  void ShowTitle(const char *s = 0);
  void ShowTemp(double T, char *s = 0);
  int  waitformouse();
  void ProcessEvents();
  void Draw();
  int Events();

  int info, topview, fh, show_connections;
  void Info(double T, double L, int antal, double acctry);
  void Info(int rader, char *s);

  double xmin, xmax, ymin, ymax;

  double xp, yp, zp;		// View...~
  double xm, ym;		// Max x,y to fit on window.
  double xmm,ymm;
  double lll;			// System size.

  void setview(double xmin, double ymin, double xmax, double ymax);

  void setview(double L);
  void tr3D(double x, double y, double z, double &xw, double &yw);
  void tr3Dp(double x, double y, double z, int &xi, int &yi);
  void dr3Dbox(double L, double Z);
  void plot3D(double x1,double y1,double z1,double x2,double y2,double z2,
	      GC gc3 = 0);
  void plot3Dsz(double x1,double y1,double z1,double x2,double y2,double z2,
		double size, GC gc3 = 0);
  //  void Redraw();
};

#endif
