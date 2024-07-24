/* 
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_CSG_ML_H
#define _VOTCA_CSG_ML_H
#pragma once

// Boost includes
#include <boost/serialization/version.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

// VOTCA includes
#include <votca/tools/eigen.h>

// Local VOTCA includes
#include "votca/csg/nblistgrid.h"
#include "votca/csg/beadlist.h"

namespace votca { namespace csg {

/**
   \brief Machine Learning Class

   Implements a machine learning class to do kernel based regression for force or energy matching.
   It stores all data needed to either

 */

class ML
{
public:

    virtual ~ML();

    /**
     * \brief constructor without abjusting the size of the objects
     * 
    */   
    ML();
    
    /**
     * \brief constructor
     * 
     * \param bead_number is the total number of beads for this ML parametrization
     * \param struct_number is the total number of structural units of the ML model
     * \param sigma is the hyperparameter sigma
     * \param threebody is the threebody flag
     * 
    */
    ML(int bead_number, int struct_number, Eigen::Vector3d sigma, bool threebody);//, double lambda);

     /**
     * \brief Resizes the vectors and matrices, if parameters have changed. Existing entries are preserved.
     * 
     * \param bead_number is the total number of beads for this ML parametrization
     * \param struct_number is the total number of structural units of the ML model
     * 
    */
    void Resize(int bead_number, int struct_number);

    /**
     * \brief Evaluates the matrix L*K to evaluate the coefficients
     * 
     * returns L*K (mapping matrix times kernel)
    */
    const Eigen::MatrixXd Evaluate_L_K_L_T();

    /**
     * \brief Predicts values for structural units according to trained ML model
     * 
     * \param descriptor_pred is vector of descriptors (structures) for the ML prediction
     * 
     * \param pred vector with the predicted quantities, as e.g. energy, according to the descriptors
    */
    void PredictStruct(Eigen::VectorXd &descriptor_pred, Eigen::VectorXd &pred,Eigen::MatrixXd &L_K_pred);

    /**
     * \brief Predicts values for structural units according to trained ML model
     * 
     * \param descriptor_pred is vector of descriptors (structures) for the ML prediction
     * 
     * \param pred vector with the predicted quantities, as e.g. energy, according to the descriptors
    */
    void PredictStruct_Threebody(Eigen::VectorXd &descriptor_pred, Eigen::VectorXi &descriptor_pred_number, Eigen::VectorXd &pred,Eigen::MatrixXd &L_K_pred);

    /**
     * \brief Predicts values for structural units according to trained ML model
     * 
     * \param descriptor_pred is vector of descriptors (structures) for the ML prediction
     * 
     * \param pred vector with the predicted quantities, as e.g. energy, according to the descriptors
    */
    void PredictStruct(Eigen::VectorXd &descriptor_pred, Eigen::VectorXd &pred);

    /**
     * \brief Predicts values for structural units according to trained ML model
     * 
     * \param descriptor_pred is vector of descriptors (structures) for the ML prediction
     * 
     * \param pred vector with the predicted quantities, as e.g. energy, according to the descriptors
    */
    void PredictStruct_Threebody(Eigen::VectorXd &descriptor_pred, Eigen::VectorXi &descriptor_pred_number, Eigen::VectorXd &pred);

    /**
     * \brief Predicts values for beads according to trained ML model
     * 
     * \param descriptor_pred is vector of descriptors (structures) for the ML prediction
     * 
     * \param mapping_pred mapping matrix for the test configuration to do the prediction
     * \param pred vector with the predicted quantities per bead, as e.g. energy
    */
    void PredictBead(Eigen::VectorXd &descriptor_pred, Eigen::MatrixXd &mapping_pred, Eigen::VectorXd &pred_beads);

    /**
     * \brief Predicts values for beads according to trained ML model
     * 
     * \param descriptor_pred is vector of descriptors (structures) for the ML prediction
     * 
     * \param mapping_pred mapping matrix for the test configuration to do the prediction
     * \param pred vector with the predicted quantities per bead, as e.g. energy
    */
    void PredictBead_Threebody(Eigen::VectorXd &descriptor_pred, Eigen::VectorXi &descriptor_pred_number, Eigen::MatrixXd &mapping_pred, Eigen::VectorXd &pred_beads);

    /**
     * get the threebody flag
     *
     * \return threebody
     */
    const bool getThreebody() const { return _threebody; }

    /**
     * set the threebody flag
     * 
     * \param threebody
     */
    void setThreebody(bool threebody) { _threebody=threebody; }
    
    /**
     * get the hyperparameter sigma
     *
     * \return sigma
     */
    const Eigen::Vector3d getSigma() const { return _sigma; }

    /**
     * set the hyperparameter sigma
     * 
     * \param sigma
     */
    void setSigma(Eigen::Vector3d sigma) { _sigma=sigma; }

    /**
     * get the current bead number
     *
     * \return _bead_number
     */
    const int getBeadNumber() const { return _bead_number; }

    /**
     * get the current bead number
     *
     * \return _struct_number
     */
    const int getStructNumber() const { return _struct_number; }

    /**
     * \brief sets value at the specific position of the mapping matrix, if not specified, the weighting is 1.0
     * 
     * \param positions in mapping matrix
     * \param scale of the weighting, if not specified = 1.0
    */
    void setMappingMatrix(int position1, int position2, double scale=1.0) { _L(position1,position2)=scale; }

    /**
     * \brief adds value to the specific position of the mapping matrix, if not specified, the weighting is 1.0
     * 
     * \param positions in mapping matrix
     * \param scale of the weighting, if not specified = 1.0
    */
    void addMappingMatrix(int position1, int position2, double scale=1.0) { _L(position1,position2)+=scale; }

    /**
     * \brief gets value at the specific position of the mapping matrix
     * 
     * \return _L(position1,position2)
     * 
    */
    const double getMappingMatrix(int position1, int position2) {return _L(position1,position2); }

    /**
     * \brief returns the whole mapping matrix
     * 
     * \return _L(position1,position2)
     * 
    */
    const Eigen::MatrixXd getMappingMatrix() { return _L; }

    /**
     * \brief sets value to descriptor
     * 
     * \param position in descriptor
     * \param val to be added
    */
    void setDescriptor(int position, double val) { _descriptor(position)=val; }

    /**
     * get the value of the descriptor at the position (position)
     *
     * \return _descriptor(position)
     */
    const double getDescriptor(int position) { return _descriptor(position); }

    /**
    */
    void setDescriptorNumber(int position, int val) { _descriptor_number(position)=val; }

    /**
     */
    const int getDescriptorNumber(int position) { return _descriptor_number(position); }

    /**
     * \brief sets value to vector of coefficients alpha
     * 
     * \param position in alpha
     * \param val to be added
    */
    void setCoefficient(int position, double val) { _alpha(position)=val; }

    /**
     * get the value of the vector of coefficients alpha at the position (position)
     *
     * \return _alpha(position)
     */
    const double getCoefficient(int position) { return _alpha(position); }

    //void EvaluateTestKernel(Eigen::VectorXd &descriptor_pred,Eigen::MatrixXd &K_pred);

private:

    /// \brief number of structural units in the system, e.g. pairs, triplets, ...
    bool _threebody;
    
    /// \brief number of structural units in the system, e.g. pairs, triplets, ...
    int _struct_number;

    /// \brief number of CG beads. The mapping matrix _L connects between the beads and the structural units.
    int _bead_number;

    /// \brief number of nonzero entries in mapping matrix _L, needed for serialization
    int _container_number;

    /// \brief number of entries in descriptor, needed for serialization
    int _descriptor_size;

    /// \brief hyperparameters sigma
    Eigen::Vector3d _sigma;
    Eigen::MatrixXi _container;
    //matrix to store value of mapping matrix != 1.0
    Eigen::VectorXd _container_value;

    /// \brief The Mapping Matrix _L to map between beads and descriptors
    Eigen::MatrixXd _L;
    /// \brief Vector to store descriptors
    Eigen::VectorXd _descriptor;
    /// \brief Vector to store which type the atom is
    Eigen::VectorXi _descriptor_number; 
    /// \brief Vector to store coefficients
    Eigen::VectorXd _alpha;

    /**
     * \brief Evaluates the Kernel Matrix, given the descriptor and hyperparameters
    */    
    void EvaluateKernel(Eigen::VectorXd &descriptor1, Eigen::VectorXd &descriptor2, Eigen::MatrixXd &K);

    /**
     * \brief Evaluates the Kernel Matrix, given the descriptor and hyperparameters
    */    
    void EvaluateKernel_Threebody(Eigen::VectorXd &descriptor1, Eigen::VectorXi &descriptor_number1, Eigen::VectorXd &descriptor2, Eigen::VectorXi &descriptor_number2, Eigen::MatrixXd &K);
    
    // Allow serialization to access non-public data members
    friend class boost::serialization::access;

    // serialization itself (template implementation stays in the header)
    template<typename Archive>
    void serialize(Archive& ar, const unsigned int version) {
        //serialize and deserialize values
        ar & _struct_number;
        ar & _bead_number;
        ar & _descriptor_size;
        for (int i=0; i<3; ++i){
            ar & _sigma(i);
        }

        //in case of deserialisation resize Eigen vector for coefficients
        if (_alpha.rows() != _bead_number*3){
            _alpha.conservativeResize(_bead_number*3);
        }
        //in case of deserialisation resize Eigen vector for descriptors
        if (_descriptor.rows() != _descriptor_size){
            _descriptor.conservativeResize(_descriptor_size);
        }
        //in case of deserialisation resize Eigen vector for descriptors
        if (_descriptor_number.rows() != _struct_number){
            _descriptor_number.conservativeResize(_struct_number);
        }
        //case of deserialisation
        if ( (_L.rows() != _bead_number*3) && (_L.cols() != _struct_number*3) ) {
            //first get
            ar & _container_number;
            //then resize mapping matrix
            _L=Eigen::MatrixXd::Zero(_bead_number*3,_struct_number*3);
            _container=Eigen::MatrixXi::Zero(_container_number,2);
            _container_value=Eigen::VectorXd::Zero(_container_number);
            for (int i = 0; i < _container_number; ++i){
                ar & _container(i,0);
                ar & _container(i,1);
                ar & _container_value(i);
                //_L(_container(i,0),_container(i,1))=1.0;
                _L(_container(i,0),_container(i,1))=_container_value(i);
            }
        }
        //case of serialization
        else{
            int count = 0;
            for (int i = 0; i < _bead_number*3; ++i){
                for (int j = 0; j < _struct_number*3; ++j){
                    if (_L(i,j)!=0)
                        ++count;
                }
            }
            _container_number=count;
            ar & _container_number;
            _container=Eigen::MatrixXi::Zero(_container_number,2);
            _container_value=Eigen::VectorXd::Zero(_container_number);
            count=0;
            for (int i = 0; i < _bead_number*3; ++i){
                for (int j = 0; j < _struct_number*3; ++j){
                    //check if there is an entry in the mapping matrix
                    if (_L(i,j)!=0) {
                        _container(count,0)=i;
                        _container(count,1)=j;
                        _container_value(count)=_L(i,j);
                        ++count;
                    }
                }
            }
            for (int i = 0; i < _container_number; ++i){
                ar & _container(i,0);
                ar & _container(i,1);
                ar & _container_value(i);
            }
        }

        for (int i = 0; i < _descriptor_size; ++i){
            ar & _descriptor(i);
        }

        for (int i = 0; i < _struct_number; ++i){
            ar & _descriptor_number(i);
        }

        for (int i = 0; i < _bead_number*3; ++i){
            ar & _alpha(i);
        }
    }// end of serialization
};    

}}

#endif /* _VOTCA_CSG_ML_H */

