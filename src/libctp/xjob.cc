#include <votca/ctp/xjob.h>

namespace votca { namespace ctp {
    

vector<XJob*> XJOBS_FROM_TABLE(const string &job_file, Topology *top) {
    
    vector<XJob*> xjobs;
    
    QMNBList &nblist = top->NBList();

    std::string line;
    std::ifstream intt;
    intt.open(job_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

// Sample line
// # JOB_ID TAG  PAIR_ID SEG1_ID SEG1_NAME SEG1_MPS SEG2_ID SEG2_NAME SEG2_MPS  (TYPE  SITE)
//   1      E_CT 3819    182     C60       c60.mps  392     DCV       dcv.mps   (site  392)

            int jobId       = boost::lexical_cast<int>(split[0]);
            string tag      = split[1];
            int pairId      = boost::lexical_cast<int>(split[2]);

            int seg1Id      = boost::lexical_cast<int>(split[3]);
            string seg1Name = split[4];
            string seg1mps  = split[5];

            int seg2Id      = boost::lexical_cast<int>(split[6]);
            string seg2Name = split[7];
            string seg2mps  = split[8];

            string job_type = "pair";
            int    energy_site_id = -1;
            if (split.size() == 11) {
                job_type = split[9];
                energy_site_id = boost::lexical_cast<int>(split[10]);
            }

            Segment *seg1   = top->getSegment(seg1Id);
            Segment *seg2   = top->getSegment(seg2Id);

            xjobs.push_back(new XJob(jobId,  tag,    job_type, energy_site_id,
                                      pairId, seg1Id, seg2Id,   seg1mps,
                                      seg2mps, top));

            xjobs.back()->setType(job_type, energy_site_id);

        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << job_file << endl;
           throw runtime_error("Please supply input file.");           }

    cout << endl 
         << "... ... ... Registered " << xjobs.size() << " jobs. "
         << flush;
    
    return xjobs;    
}

    
}}