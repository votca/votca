/* 
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#ifndef VOTCA_XTP_JOB_H
#define	VOTCA_XTP_JOB_H

#include <iostream>
#include <fstream>
#include <string>
#include <votca/tools/property.h>

using namespace std;
using namespace votca::tools;

namespace votca { namespace xtp {
    

    
    
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
    Job(int id, string &tag, string &input, string stat);
    Job(int id, string &tag, Property &input, JobStatus stat);
   ~Job() {;}
    
    string    ConvertStatus(JobStatus) const;
    JobStatus ConvertStatus(string) const;
   
    class JobResult
    {
    public:
        
        JobResult() { ; }
        
        void setStatus(JobStatus stat) { _status = stat; }
        void setStatus(string stat) { assert(false); }
        void setOutput(string output) 
            { _has_output = true; _output = Property().add("output", output); }
        void setOutput(Property &output)
            { _has_output = true; _output = output.get("output"); }
        void setError(string error) { _has_error = true; _error = error; }
        
        JobStatus _status;
        Property _output;
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
    Property &getInput() { return _input; }    
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
    void setOutput(string output) { _output = Property().add("output", output); _has_output = true; }
   
    const string &getHost() const { assert(_has_host); return _host; }
    const string &getTime() const { assert(_has_time); return _time; }
    const Property &getOutput() const { assert(_has_output); return _output; }
    const string &getError() const { assert(_has_error); return _error; }

protected:

    
     // Defined by user
     int _id;
     string _tag;    
     JobStatus _status;
     int _attemptsCount;
     Property _input;
     string _sqlcmd;
    
     // Generated during runtime
     string _host;
     bool   _has_host;
     string _time;
     bool   _has_time;
     Property _output;
     bool   _has_error;
     bool   _has_output;
     string _error;
     bool   _has_sqlcmd;
};




    
    
    
    
    
}}


#endif
