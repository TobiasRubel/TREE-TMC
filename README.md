TREE-QMC
--------

To build TREE-QMC, use commands:
```
git clone https://github.com/yhhan19/TREE-QMC.git
cd TREE-QMC/MQLib
make
cd ..
 g++ -std=c++11 -O2 -I MQLib/include -o TREE-QMC src/TREE-QMC.cpp MQLib/bin/MQLib.a
```

To run TREE-QMC, use command:
```
./TREE-QMC -i <input file> -o <output name>
```

To see TREE-QMC usage options, use command:
```
./TREE-QMC -h
```
