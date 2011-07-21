#include "qmdatabase.h"

void QMDatabase::onCreate()
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
      
    // table for conjugated segments
    Exec("CREATE TABLE conjsegs ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "name TEXT NOT NULL,"
        "type TEXT NOT NULL,"
        "molecule INT NOT NULL,"
        "frame INT NOT NULL," // TODO: this is bad
        "occ REAL)"); 

    // additional properties of conjugated segments
    Exec("CREATE TABLE conjseg_properties ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "conjseg INTEGER NOT NULL,"
        "key TEXT NOT NULL,"
        "value REAL NOT NULL)");

    // table for rigid fragments
    Exec("CREATE TABLE rigidfrags ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "frame INT NOT NULL,"
        "name TEXT NOT NULL,"
        "symmetry INT NOT NULL,"
        "type TEXT NOT NULL,"
        "resnr INT NOT NULL,"
        "mass REAL NOT NULL,"
        "charge REAL NOT NULL,"
        "molecule INT NOT NULL,"
        "conjseg_id INT NOT NULL,"
        "conjseg_index INT NOT NULL,"
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
    Exec("CREATE TABLE rigidfrag_properties ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "fragmenttid INTEGER NOT NULL,"
        "key TEXT NOT NULL,"
        "value REAL NOT NULL)");
               
    // table for pairs
    Exec("CREATE TABLE pairs ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "conjseg1 INT NOT NULL,"
        "conjseg2 INT NOT NULL,"
        "rate12 REAL NOT NULL,"
        "rate21 REAL NOT NULL,"
        "r_x REAL NOT NULL,"
        "r_y REAL NOT NULL,"
        "r_z REAL NOT NULL,"
        "deleted INT DEFAULT 0)");

    // additional properties of pairs
    Exec("CREATE TABLE pair_properties ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "pair INTEGER NOT NULL,"
        "key TEXT NOT NULL,"
        "value REAL NOT NULL)");

    Exec("CREATE TABLE pair_integrals ("
        "pair INTEGER NOT NULL,"
        "num INTEGER NOT NULL,"
        "J REAL NOT NULL)");


/*
    // delete molecules if frame is deleted
    Exec("CREATE TRIGGER trig_delete_frame BEFOR DELETE ON frames "
            "FOR EACH ROW BEGIN "
            "DELETE FROM molecules WHERE molecules.frame = OLD._id;"
            " END");

    // delete segment properties + fragments + pairs if segment is deleted
    Exec("CREATE TRIGGER trig_delete_conjseg BEFOR DELETE ON conjseg "
            "FOR EACH ROW BEGIN "
            "DELETE FROM conjseg_properties WHERE conjseg_properties.conjsegment = OLD._id;"
            "DELETE FROM rigidfrag WHERE conjseg_id= OLD._id;"
            "DELETE FROM pairs WHERE pairs.conjseg1 = OLD._id OR pairs.conjseg2 = OLD._id;"
            " END");

    // delete fragment properties if fragment is deleted
    Exec("CREATE TRIGGER trig_delete_fragment BEFOR DELETE ON rigidfrag "
            "FOR EACH ROW BEGIN "
            "DELETE FROM rigidfrag_properties WHERE fragment_properties.fragment = OLD._id;"
            "DELETE FROM conjseg WHERE segment_id=segmentid.frame = OLD._id;"
            " END");

    // delete pair property if property is deleted
    Exec("CREATE TRIGGER trig_delete_pair BEFOR DELETE ON pairs "
            "FOR EACH ROW BEGIN "
            "DELETE FROM pair_properties WHERE pair_properties.pair = OLD._id;"
            " END");

*/
}
