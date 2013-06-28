#ifndef VOTCA_CTP_JOB_H
#define	VOTCA_CTP_JOB_H

#include <iostream>
#include <fstream>
#include <string>
#include <votca/tools/property.h>

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {
    

    
    
class Job
{

public:
    
    enum JobStatus 
    { 
        AVAILABLE, 
        ASSIGNED, 
        FAILED, 
        COMPLETE
    };
    

    Job(Property *prop);
    Job(int id, string &tag, string &input, string &stat);
   ~Job() {;}
    
    string    ConvertStatus(JobStatus) const;
    JobStatus ConvertStatus(string) const;
   
    class JobResult
    {
    public:
        
        JobResult() { ; }
        
        void setStatus(JobStatus stat) { _status = stat; }
        void setStatus(string stat) { assert(false); }
        void setOutput(string output) { _has_output = true; _output = output; }
        void setError(string error) { _has_error = true; _error = error; }
        
        JobStatus _status;
        string _output;
        bool _has_output;
        string _error;
        bool _has_error;
    };

    void Reset();
    void ToStream(ofstream &ofs, string fileformat);
    void UpdateFrom(Job *ext);
    void SaveResults(JobResult *res);
   
    int getId() const { return _id; }
    string getTag() const { return _tag; }
    string getInput() const { return _input; }    
    const JobStatus &getStatus() const { return _status; }
    string getStatusStr() const { return ConvertStatus(_status); }
    
    bool hasHost() const { return _has_host; }
    bool hasTime() const { return _has_time; }
    bool hasOutput() const { return _has_output; }
    bool hasError() const { return _has_error; }
    
    bool isAvailable() const { return (_status == AVAILABLE) ? true : false; }
    bool isAssigned() const { return (_status == ASSIGNED) ? true : false; }
    bool isFailed() const { return (_status == FAILED) ? true : false; }
    bool isComplete() const { return (_status == COMPLETE) ? true : false; }
    bool isFresh() const { return (_attemptsCount < 1) ? true : false; }
    
    void setStatus(JobStatus stat) { _status = stat; }
    void setStatus(string stat) { _status = ConvertStatus(stat); }
    void setTime(string time) { _time = time; _has_time = true; }
    void setHost(string host) { _host = host; _has_host = true; }
    void setOutput(string output) { _output = output; _has_output = true; }
   
    const string &getHost() const { assert(_has_host); return _host; }
    const string &getTime() const { assert(_has_time); return _time; }
    const string &getOutput() const { assert(_has_output); return _output; }
    const string &getError() const { assert(_has_error); return _error; }

protected:

    
     // Defined by user
     int _id;
     string _tag;    
     JobStatus _status;
     int _attemptsCount;
     string _input;
     string _sqlcmd;
     bool   _has_sqlcmd;
    
     // Generated during runtime
     string _host;
     bool   _has_host;
     string _time;
     bool   _has_time;
     string _output;
     bool   _has_output;
     string _error;
     bool   _has_error;
};




    
    
    
    
    
}}


#endif