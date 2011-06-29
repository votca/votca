#include "qmdatabase.h"

QMDatabase::onCreate()
{
   // table for frames
    Exec(
        "CREATE TABLE frames ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "time REAL NOT NULL,"
        "step INT NOT NULL,"
        "box11 REAL NOT NULL,"
        "box12 REAL NOT NULL,"
        "box13 REAL NOT NULL,"
        "box21 REAL NOT NULL,"
        "box22 REAL NOT NULL,"
        "box23 REAL NOT NULL,"
        "box31 REAL NOT NULL,"
        "box32 REAL NOT NULL,"
        "box33 REAL NOT NULL)");
   
    // table for molecules
    Exec(
        "CREATE TABLE molecules ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL)");
       
    // table for conjugated segments
    Exec(
        "CREATE TABLE crgunits ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL,"
        "energy REAL NOT NULL,"
        "occ REAL NOT NULL)");
           
    // table for rigid fragments
    Exec(
        "CREATE TABLE beads ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL,"
        "symmetry INT NOT NULL,"
        "type TEXT NOT NULL,"
        "resnr INT NOT NULL,"
        "mass REAL NOT NULL,"
        "charge REAL NOT NULL,"
        "crgunit INT NOT NULL,"
        "crgunit_index INT NOT NULL,"
        "pos_x REAL NOT NULL,"
        "pos_y REAL NOT NULL,"
        "pos_z REAL NOT NULL,"
        "u_x REAL NOT NULL,"
        "u_y REAL NOT NULL,"
        "u_z REAL NOT NULL,"
        "v_x REAL NOT NULL,"
        "v_y REAL NOT NULL,"
        "v_z REAL NOT NULL,"
        "molid INT NOT NULL)");
               
    // table for pairs
    Exec(
        "CREATE TABLE pairs ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "crgunit1 INT NOT NULL,"
        "crgunit2 INT NOT NULL,"
        "rate12 REAL NOT NULL,"
        "rate21 REAL NOT NULL,"
        "r_x REAL NOT NULL,"
        "r_y REAL NOT NULL,"
        "r_z REAL NOT NULL)");
}
