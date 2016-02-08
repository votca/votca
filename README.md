This is VOTCA's next generation buildutil. 

Usage:
```
git clone --recursive -b global https://github.com/votca/buildutil.git
cd buildutil
mkdir build
cd build
cmake -DBUILD_CSGAPPS=ON -DBUILD_CTP=ON -DBUILD_XTP=ON ..
make -j5
```

Further information on VOTCA can be found at http://www.votca.org
