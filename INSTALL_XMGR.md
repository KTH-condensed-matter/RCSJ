# Installing XMGR

This repository was written to use xmgr, which is currently no longer maintained and therefore somewhat difficult to install. While we are migrating this repository to use something more modern this guide can be used to install xmgr from source.

This instruction builds on the Github repository [mlund/xmgr-resurrection](https://github.com/mlund/xmgr-resurrection). Please refer to it if support is needed.

## Installation

Install library dependencies. On Debian/Ubuntu run

```shell script
sudo apt install libice-dev libx11-dev libmotif-dev libxmu-dev libxpm-dev cmake
```

Clone the repository

```shell script
git clone git@github.com:mlund/xmgr-resurrection.git
```

Configure build

```shell script
cmake . -DCMAKE_INSTALL_PREFIX=/usr/local
```

Build project

```shell script
make
```

Install xmgr

```shell script
sudo make install
```

To use xmgr you need to add it to path. This can be done by adding the following line in the `.profile` or `.bash_profile` (or similar depending on your shell)

```shell script
PATH="$PATH:/usr/local/xmgr/bin"
```

Restart the shell or source the file updated. The command xmgr is now available in terminal.
