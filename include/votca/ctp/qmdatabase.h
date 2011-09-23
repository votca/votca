#ifndef __VOTCA_CTP_QMDATABASE_H
#define	__VOTCA_CTP_QMDATABASE_H

#include <votca/tools/database.h>

using namespace votca::tools;

namespace votca { namespace ctp {

/**
 * \brief the state database
 *
 * This class contains management of state databases. It creates
 * new databases and opens or upgrades existing ones.
 */
class QMDatabase :
    public Database
{
public:
    
    /**
     * \brief Create the database scheme
     *
     * This function is called when a database is created.
     * All tables and triggers should be added here.
     */
    void onCreate();
};

}}

#endif	/* QMDATABASE_H */

