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

#ifndef _CSG_ML_H
#define	_CSG_ML_H

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/csg/csgapplication.h"
#include "votca/csg/trajectoryreader.h"
#include "votca/csg/ml.h"

using namespace votca::csg;

using namespace std;

/**
    \brief machine learning algorithm using kernel based regression
 **/

class CGMachineLearning
    : public CsgApplication
{
public:

    string ProgramName() override { return "csg_ml"; }
    void HelpText(ostream &out) {
        out << "Performs kernel-based machine learning (ML) for non-bonded two-body and three-body interactions. Can be used either for training or testing of a ML model, or to evaluate a three-body force table.";
    }

    bool DoTrajectory() override {return true;}
    bool DoMapping() override {return true;}

    void Initialize(void) override;
    bool EvaluateOptions() override;

    /// \brief called before the first frame
    void BeginEvaluate(Topology *top, Topology *top_atom) override;
    /// \brief called after the last frame
    void EndEvaluate() override;
    /// \brief called for each frame which is mapped
    void EvalConfiguration(Topology *conf, Topology *conf_atom = nullptr) override;
    /// \brief load options from the input file
    void LoadOptions(const string &file);

protected:
    /// \brief structure, which contains ML object with related parameters, so far only for nonbonded interactions
    struct MLInfo {
        /// \brief constructor
        MLInfo(int index, votca::tools::Property *options);

        /// \brief interaction index
        int MLIndex;
        /// \brief ML Name
        string MLName;        
        /// \brief types of beads involved
        string type1, type2, type3;

        /// \brief pointer to Property object to handle input options
        votca::tools::Property *_options;

        /// \brief ML object
        ML *MLObject;

        /// \brief true if non-bonded interaction is threebody interaction
        bool threebody;
        /// \brief true if non-bonded threebody interaction has symmetric tables for binning and output (type2 == type3)
        bool threebody_symmetric;
        /// \brief minimum value for force output on sample triplets
        double min_out;
        /// \brief maximum value for force output on sample triplets
        double max_out;
        /// \brief dx for force output on sample triplets
        double dx_out;
        /// \brief number of points for force output on sample triplets
        int num_out;
	/// \brief range for applying cutoff function
        double d;
        /// \brief accumulated test error
        double test_error_sum;

        Eigen::MatrixXd result;
        Eigen::MatrixXd resulttheta;
        /// \brief accuracy of the final result
        Eigen::MatrixXd error;
        Eigen::MatrixXd errortheta;

        /// \brief yes if output table should be evaluated
        bool output_table;
        /// \brief dx for force output in output table
        double dx_table;
        /// \brief number of points for force output in output table
        int num_table;

        /// \brief vectors to store output values for force table
        Eigen::VectorXd output1_1;
        Eigen::VectorXd output1_2;
        Eigen::VectorXd output2_1;
        Eigen::VectorXd output2_2;
        Eigen::VectorXd output3_1;
        Eigen::VectorXd output3_2;

        //Extra variables for binning

        /// \brief true if binning is used for this interaction
        bool binning;
        /// \brief minimum value for pair distances of second grid
        double min_out2;
        /// \brief minimum theta value for angular grid        
        double min_theta;
	/// \brief number of bins
	int num_bins;
	/// \brief number of bins
	int num_theta;
	/// \brief number of bins for second grid
	int num_bins2;
	/// \brief dx for bins
	double dx_bins;
	/// \brief dx for bins of second grid
	double dx_bins2;
	/// \brief scaling for smearing width when assigning to bins
	double smear_scale;  
        /// \brief determines if the binning should be done with Gaussian smearing.
        bool gaussian_smearing;

        //Extra vector for storing the binned pairs or triplets
        Eigen::VectorXd binned_structures;
    };
    /// \brief Property object to handle input options
    votca::tools::Property _options;

    /// \brief list of non-bonded interactions
    std::vector<votca::tools::Property *> _nonbonded;

    using MLContainer = vector<MLInfo>;
    /// \brief vector of SplineInfo * for all interactions
    MLContainer _mls;

    /// \brief matrix used to store matrix for ML equations
    Eigen::MatrixXd _L_K_L_T;
    /// \brief vector used to store reference forces or energies on CG beads (from atomistic simulations)
    Eigen::VectorXd _b;
    /// \brief vector used to store mapping between randomly chosen bead indices and rows in mapping matrix
    Eigen::VectorXi _bead_to_row;
    /// \brief vector used to temporarily store the coefficients
    Eigen::VectorXd _x;

    /// \brief Counter for trajectory frames
    int _frame_counter;
    /// \brief Number of CG beads
    int _nbeads;
    /// \brief Number of frames used to train or test
    int _nframes;
    /// \brief Number of beads for frames used 
    int _nbeads_per_frame;

    /// \brief accuracy for evaluating the difference in bead positions
    double _dist;

    /// \brief regularization parameter for linear equations
    double _lambda;

    /// \brief determines if the code is run in training or testing mode.
    bool _train;

    /// \brief determines if the beads of one frame are randomly selected or in order of their numbering.
    bool _random_selection;

    bool has_existing_forces_;

    /// \brief Solve ML equations and accumulate the result
    void AccumulateDataTrain();
    // combine test results for different frames
    void AccumulateDataTest();

    /// \brief For each trajectory frame sets mapping, descriptors and matrices
    void EvalNonbondedTrain(Topology *conf, MLInfo *mlinfo);
    /// \brief For each trajectory frame sets mapping, descriptors and matrices for threebody interactions
    void EvalNonbondedTrain_Threebody(Topology *conf, MLInfo *mlinfo);
    /// \brief For each trajectory frame does testing
    void EvalNonbondedTest(Topology *conf, MLInfo *mlinfo);
    /// \brief For each trajectory frame does testing for threebody interactions
    void EvalNonbondedTest_Threebody(Topology *conf, MLInfo *mlinfo);
    /// \brief Write results to output files
    void SerializeMLObjects();
    /// \brief Write results to output files
    void deSerializeMLObjects();
    /// \brief Write results to output files
    void WriteOutFilesTrain();

    /// \brief Evaluate triplet configurations to write to output table
    void EvaluateTable_Threebody(MLInfo *mlinfo);
    /// \brief write output table
    void WriteTable_Threebody(MLInfo *mlinfo);

    const double Calculate_fcut(double r, double r_cut, double d);

    const Eigen::Matrix3d get_Rxz(Eigen::Vector3d &vec);

    const Eigen::Matrix3d get_Rz(Eigen::Vector3d &vec);

    const Eigen::Matrix3d get_Rvx(Eigen::Vector3d &vec);

    const Eigen::Matrix3d get_Rvz(Eigen::Vector3d &vec);

    Topology top_force_;
    std::unique_ptr<TrajectoryReader> trjreader_force_;
};

#endif	/* _CSG_ML_H */
