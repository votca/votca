This is VOTCA's next generation build system.

Usage:

```
prefix=WHERE/TO/INSTALL/VOTCA
git clone --recursive https://github.com/votca/votca.git
cd votca
mkdir build
cd build
cmake -DBUILD_CSGAPPS=ON -DBUILD_CTP=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} ..
make -j5
```

Using this code via docker:

```
docker pull votca/votca
docker run -it votca/votca /bin/bash
```

For further details see:
1. [Installation](share/doc/INSTALL.md) 
2. [Further Information](http://www.votca.org)
3. [Developers Guide](share/doc/DEVELOPERS_GUIDE.md)

