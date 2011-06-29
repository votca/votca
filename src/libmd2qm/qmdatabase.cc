#include "qmdatabase.h"

QMDatabase::onCreate()
{
   // table for frames
    Exec("CREATE TABLE frames ("
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
    Exec("CREATE TABLE molecules ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL,"
        "frame INT NOT NULL)");

    // delete molecules if frame is deleted
    Exec("CREATE TRIGGER trig_delete_frame BEFOR DELETE ON frames "
            "FOR EACH ROW BEGIN "
            "DELETE FROM molecules WHERE molecules.frame = OLD._id;"
            " END");
       
    // table for conjugated segments
    Exec("CREATE TABLE segments ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL,"
        "occ REAL NOT NULL)");

    // additional properties of conjugated segments
    Exec("CREATE TABLE segment_properties ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "segmentid INTEGER NOT NULL,"
        "key TEXT NOT NULL,"
        "value REAL NOT NULL)");

        // delete molecules if frame is deleted
    Exec("CREATE TRIGGER trig_delete_frame BEFOR DELETE ON frames "
            "FOR EACH ROW BEGIN "
            "DELETE FROM molecules WHERE molecules.frame = OLD._id;"
            " END");

    // table for rigid fragments
    Exec("CREATE TABLE fragments ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL,"
        "symmetry INT NOT NULL,"
        "type TEXT NOT NULL,"
        "resnr INT NOT NULL,"
        "mass REAL NOT NULL,"
        "charge REAL NOT NULL,"
        "molid INT NOT NULL,"
        "segment_id INT NOT NULL,"
        "segment_index INT NOT NULL,"
        "pos_x REAL NOT NULL,"
        "pos_y REAL NOT NULL,"
        "pos_z REAL NOT NULL,"
        "u_x REAL NOT NULL,"
        "u_y REAL NOT NULL,"
        "u_z REAL NOT NULL,"
        "v_x REAL NOT NULL,"
        "v_y REAL NOT NULL,"
        "v_z REAL NOT NULL)");

    // additional properties of rigid fragments
    Exec("CREATE TABLE fragment_properties ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "fragmenttid INTEGER NOT NULL,"
        "key TEXT NOT NULL,"
        "value REAL NOT NULL)");
               
    // table for pairs
    Exec("CREATE TABLE pairs ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "segment1 INT NOT NULL,"
        "segment2 INT NOT NULL,"
        "rate12 REAL NOT NULL,"
        "rate21 REAL NOT NULL,"
        "r_x REAL NOT NULL,"
        "r_y REAL NOT NULL,"
        "r_z REAL NOT NULL)");

    // additional properties of pairs
    Exec("CREATE TABLE pair_properties ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "pairid INTEGER NOT NULL,"
        "key TEXT NOT NULL,"
        "value REAL NOT NULL)");
}
