/*
 *            Copyright 2009-2018 The VOTCA Development Team
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


#include <iostream>
#include <fstream>
#include <stdexcept>
#include <votca/xtp/votca_config.h>

#include <string>
#include <typeinfo>
#include <vector>
#include <boost/filesystem.hpp>
#include <type_traits>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/checkpointreader.h>
#include <votca/xtp/checkpoint.h>
namespace votca {
namespace xtp {

using namespace checkpoint_utils;
namespace bfs = boost::filesystem;

std::ostream& operator<<(std::ostream& s, CheckpointAccessLevel l){

    switch(l){
    case CheckpointAccessLevel::READ:
        s << "read"; break;
    case CheckpointAccessLevel::MODIFY:
        s << "modify"; break;
    case CheckpointAccessLevel::CREATE:
        s << "create"; break;
    }

    return s;
}

bool FileExists(const std::string& fileName){
    return bfs::exists(fileName);
}

CheckpointFile::CheckpointFile(std::string fN):
    CheckpointFile(fN, CheckpointAccessLevel::MODIFY) {};

CheckpointFile::CheckpointFile(std::string fN, CheckpointAccessLevel access)
    : _fileName(fN), _accessLevel(access){

    try {
        H5::Exception::dontPrint();

        switch (_accessLevel){
        case CheckpointAccessLevel::READ:
            _fileHandle = H5::H5File(_fileName, H5F_ACC_RDONLY);
            break;
        case CheckpointAccessLevel::CREATE:
            _fileHandle = H5::H5File(_fileName, H5F_ACC_TRUNC);
            break;
        case CheckpointAccessLevel::MODIFY:
            if (!FileExists(_fileName))
                _fileHandle = H5::H5File(_fileName, H5F_ACC_TRUNC);
            else
                _fileHandle = H5::H5File(_fileName, H5F_ACC_RDWR);

        }

    } catch (H5::Exception& error) {
        std::stringstream message;
        message << "Could not access file " << _fileName;
        message << " with permission to " << _accessLevel << "." << std::endl;

        throw std::runtime_error(message.str());
    }
};

std::string CheckpointFile::getFileName() { return _fileName; };

H5::H5File CheckpointFile::getHandle() { return _fileHandle; };

CheckpointWriter CheckpointFile::getWriter(const std::string _path){
    if (_accessLevel == CheckpointAccessLevel::READ){
        throw std::runtime_error("Checkpoint file opened as read only.");
    }

    try{
        return CheckpointWriter(_fileHandle.createGroup(_path), _path);
    } catch(H5::Exception& error) {
        try {
            return CheckpointWriter(_fileHandle.openGroup(_path), _path);
        } catch (H5::Exception& error) {
            std::stringstream message;
            message << "Could not create or open " << _fileName << ":" << _path
                    << std::endl;

            throw std::runtime_error(message.str());
        }
    }
};

CheckpointWriter CheckpointFile::getWriter(){
    return getWriter("/");
};

CheckpointReader CheckpointFile::getReader(const std::string _path){
    try{
        return CheckpointReader(_fileHandle.openGroup(_path), _path);
    } catch (H5::Exception& error){
        std::stringstream message;
        message << "Could not open " << _fileName << ":" << _path
                << std::endl;

        throw std::runtime_error(message.str());
    }
};

CheckpointReader CheckpointFile::getReader(){
    return getReader("/");
};

}  // namespace xtp
}  // namespace votca
