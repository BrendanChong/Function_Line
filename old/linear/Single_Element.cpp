#include <string>
#include <iostream>
#include <sstream>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Bitset.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#endif

    //------------------------------------------------------------------------------
    namespace moris
    {
        /* ------------------------------------------------------------------------ */
        // function to convert real value into string

        std::string
        moris_to_string(real tValue)
        {
            std::ostringstream streamObj;

            // Set precision
            streamObj << std::scientific;
            streamObj << std::setprecision(15);

            // Add value to stream
            streamObj << tValue;

            // Get string from output string stream
            return streamObj.str();
        }

        /* ------------------------------------------------------------------------ */
        // problem type: forward only or optimizationf
        bool tIsOptimization = true;

        /* ------------------------------------------------------------------------ */
        // Background mesh

        uint tNE = 68;
        bool tUseConsistentSens = true;
        bool tUseMultilinearInterpolation = false;

        /* ------------------------------------------------------------------------ */
        // Geometry parameters

        real tHeight = M_PI / 4.0;
        real tWidth = M_PI / 4.0;

        /* ------------------------------------------------------------------------ */
        // loading parameters

        // prescribed heat flux - volumetric
        real tHeatFlx = 100.0;

        /* ------------------------------------------------------------------------ */
        // material parameters

        std::string tDensity = moris_to_string(1.0);
        std::string tCapacity = moris_to_string(1.0);
        std::string tConductivity = moris_to_string(1.916229656078740e-01);

        /* ------------------------------------------------------------------------ */
        // parameters for weak enforcement of boundary conditions and ghost stabilization

        // Nitsche penalty
        std::string tNitschePenThermal = "100.0";

        // Flag to turn on/off ghost
        bool tUseGhost = false;

        /* ------------------------------------------------------------------------ */
        // Nonlinear solver parameters

        // maximum number of Newton steps
        int tNLA_max_iter = 5;

        // required drop of resdiual
        real tNLA_rel_res_norm_drop = 1e-12;

        // relaxation parameter
        real tNLA_relaxation_parameter = 1.0;

        /* ------------------------------------------------------------------------ */
        // File names

        std::string tName = "Single_Element";
        std::string tExoFile = tName + ".exo";
        std::string tSoFile = tName + ".so";
        std::string tHdf5File = tName + "_SEN.hdf5";
        std::string tGENOutputFile = tName + "_GEN.exo";

        /* ------------------------------------------------------------------------ */
        // Optimization parameters

        // FD in adjoint
        real tFEMFdEpsilon = 1.0e-8;

        // FD in sweep
        std::string tSweepFdEpsilon = "1.0e-6";

        // number of constraints (here: number of design criteria)
        uint tNumConstraints = 1;

        /* ------------------------------------------------------------------------ */
        // Output Config
        bool tUsePhases = false;
        std::string tSolidPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";
        std::string tVoidPhase = "HMR_dummy_n_p1,HMR_dummy_c_p1";
        std::string tInterface = "iside_b0_0_b1_1";
        std::string tAllPhases = tSolidPhase + "," + tVoidPhase + "," + tInterface;

        /* ------------------------------------------------------------------------ */
        void
        Func_HeatLoad(
            Matrix<DDRMat> &aPropMatrix,
            Vector<Matrix<DDRMat>> &aParameters,
            fem::Field_Interpolator_Manager *aFIManager)
        {
            aPropMatrix.set_size(1, 1);

            const Matrix<DDRMat> tCoord = aFIManager->get_IP_geometry_interpolator()->valx();

            aPropMatrix(0, 0) = tHeatFlx;
        }

        // Level-Set functions
        // ------------------------------------------------------------------------

        real
        Bottom_Line(
            const Matrix<DDRMat> &aCoordinates,
            const Vector<real> &aGeometryParameters)
        {
            real tTheta = aGeometryParameters(0);
            real tC = cos(tTheta);
            real tS = sin(tTheta);
            real tXNormal = tS;
            real tYNormal = -tC;

            real tXN = 0.0;
            real tYN = -tHeight;
            real tXCenter = tC * tXN - tS * tYN;
            real tYCenter = tS * tXN + tC * tYN;

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter);
        }

        real
        Left_Line(
            const Matrix<DDRMat> &aCoordinates,
            const Vector<real> &aGeometryParameters)
        {
            real tTheta = aGeometryParameters(0);
            real tC = cos(tTheta);
            real tS = sin(tTheta);
            real tXNormal = -tC;
            real tYNormal = -tS;

            real tXN = -tWidth;
            real tYN = 0.0;
            real tXCenter = tC * tXN - tS * tYN;
            real tYCenter = tS * tXN + tC * tYN;

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter);
        }

        real
        Right_Line(
            const Matrix<DDRMat> &aCoordinates,
            const Vector<real> &aGeometryParameters)
        {
            real tTheta = aGeometryParameters(0);
            real tC = cos(tTheta);
            real tS = sin(tTheta);
            real tXNormal = tC;
            real tYNormal = tS;

            real tXN = tWidth;
            real tYN = 0.0;
            real tXCenter = tC * tXN - tS * tYN;
            real tYCenter = tS * tXN + tC * tYN;

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter);
        }

        real
        Top_Line(
            const Matrix<DDRMat> &aCoordinates,
            const Vector<real> &aGeometryParameters)
        {
            real tTheta = aGeometryParameters(0);
            real tC = cos(tTheta);
            real tS = sin(tTheta);
            real tXNormal = -tS;
            real tYNormal = tC;

            real tXN = 0.0;
            real tYN = tHeight;
            real tXCenter = tC * tXN - tS * tYN;
            real tYCenter = tS * tXN + tC * tYN;
            real tCurrHeight = aGeometryParameters(1); // design variable

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter) - tCurrHeight;
        }

        // derivative of level set function wrt parameters
        void
        dLine_dHeight(
            const Matrix<DDRMat> &aCoordinates,
            const Vector<real> &aGeometryParameters,
            Matrix<DDRMat> &aFieldSensitivity)
        {
            aFieldSensitivity = {{0.0}};
        }

        // derivative of level set function wrt parameters
        void
        dTop_Line_dHeight(
            const Matrix<DDRMat> &aCoordinates,
            const Vector<real> &aGeometryParameters,
            Matrix<DDRMat> &aFieldSensitivity)
        {
            aFieldSensitivity = {{0.0, -1.0}};
        }

        void
        dTop_Line_dX(
            const Matrix<DDRMat> &aCoordinates,
            const Vector<real> &aGeometryParameters,
            Matrix<DDRMat> &aFieldSensitivity)
        {

            real tTheta = aGeometryParameters(0);
            real tC = cos(tTheta);
            real tS = sin(tTheta);
            real tXNormal = -tS;
            real tYNormal = tC;

            aFieldSensitivity(0) = tXNormal;
            aFieldSensitivity(1) = tYNormal;
        }

        /* ------------------------------------------------------------------------ */
        // Phase assignement

        uint get_phase_index(const Bitset<4> &aGeometrySigns)
        {
            // Bottom void for dirichlet boundary condition
            if (aGeometrySigns.test(0))
            {
                return 2;
            }
            // Top void for Neumann boundary condition
            if (aGeometrySigns.test(3))
            {
                return 1;
            }
            // Side voids (unused)
            if (aGeometrySigns.test(1) or aGeometrySigns.test(2))
            {
                return 3;
            }
            // Solid phase
            else
            {
                return 0;
            }
        }
        /* ------------------------------------------------------------------------ */

        Matrix<DDSMat>
        get_constraint_types()
        {
            Matrix<DDSMat> tConstraintTypes(tNumConstraints, 1, 1);

            return tConstraintTypes;
        }

        Matrix<DDRMat>
        compute_objectives(
            const Vector<real> &aADVs,
            const Vector<real> &aCriteria)
        {
            Matrix<DDRMat> tObjectives = {{aCriteria(0)}};

            return tObjectives;
        }

        /* ------------------------------------------------------------------------ */

        Matrix<DDRMat>
        compute_dobjective_dadv(
            const Vector<real> &aADVs,
            const Vector<real> &aCriteria)
        {
            Matrix<DDRMat> tDObjectiveDADV(1, aADVs.size(), 0.0);

            return tDObjectiveDADV;
        }

        /* ------------------------------------------------------------------------ */

        Matrix<DDRMat>
        compute_dobjective_dcriteria(
            const Vector<real> &aADVs,
            const Vector<real> &aCriteria)
        {
            Matrix<DDRMat> tDObjectiveDCriteria(1, aCriteria.size(), 0.0);

            tDObjectiveDCriteria(0) = 1.0;

            return tDObjectiveDCriteria;
        }

        /* ------------------------------------------------------------------------ */

        Matrix<DDRMat>
        compute_constraints(
            const Vector<real> &aADVs,
            const Vector<real> &aCriteria)
        {
            Matrix<DDRMat> tConstraints(tNumConstraints, 1);

            // Mass flux constraint
            tConstraints(0) = aCriteria(0);

            return tConstraints;
        }

        /* ------------------------------------------------------------------------ */

        Matrix<DDRMat>
        compute_dconstraint_dadv(
            const Vector<real> &aADVs,
            const Vector<real> &aCriteria)
        {
            Matrix<DDRMat> tDConstraintDADV(tNumConstraints, aADVs.size(), 0.0);

            return tDConstraintDADV;
        }

        /* ------------------------------------------------------------------------ */

        Matrix<DDRMat>
        compute_dconstraint_dcriteria(
            const Vector<real> &aADVs,
            const Vector<real> &aCriteria)
        {
            Matrix<DDRMat> tDConstraintDCriteria(tNumConstraints, aCriteria.size(), 0.0);

            tDConstraintDCriteria(0, 0) = 1.0;

            return tDConstraintDCriteria;
        }

        /* ------------------------------------------------------------------------ */
        // function to dtermine when to oupt results

        bool
        Output_Criterion(tsa::Time_Solver *aTimeSolver)
        {
            return true;
        }

        /* ------------------------------------------------------------------------ */

        void
        OPTParameterList(Module_Parameter_Lists &aParameterLists)
        {
            aParameterLists.set("is_optimization_problem", tIsOptimization);
            aParameterLists.set("problem", "user_defined");
            aParameterLists.set("library", tSoFile);

            aParameterLists(OPT::ALGORITHMS).add_parameter_list(opt::Optimization_Algorithm_Type::SWEEP);
            aParameterLists.set("hdf5_path", "SEN_" + std::to_string(tNE) + ".hdf5");
            aParameterLists.set("num_evaluations_per_adv", "30,1"); // sweep over theta (first ADV)
            aParameterLists.set("finite_difference_type", "all");
            aParameterLists.set("finite_difference_epsilons", tSweepFdEpsilon);
            aParameterLists.set("finite_difference_adv_indices", "1"); // sensitivity over height
            aParameterLists.set("evaluate_objective_gradients", true);
            aParameterLists.set("evaluate_constraint_gradients", false);
        }

        /* ------------------------------------------------------------------------ */

        void
        HMRParameterList(Module_Parameter_Lists &aParameterLists)
        {
            aParameterLists.set("number_of_elements_per_dimension", tNE, tNE);
            aParameterLists.set("domain_dimensions", 3.0, 3.0);
            aParameterLists.set("domain_offset", -1.5, -1.5);

            aParameterLists.set("refinement_buffer", 1);
            aParameterLists.set("staircase_buffer", 1);
            aParameterLists.set("pattern_initial_refinement", (uint)0, (uint)0);

            aParameterLists(HMR::LAGRANGE_MESHES).add_parameter_list();
            aParameterLists.set("order", 1);
            aParameterLists.set("pattern_index", 1);

            aParameterLists(HMR::BSPLINE_MESHES).add_parameter_list();
            aParameterLists.set("orders", 2);
            aParameterLists.set("pattern_index", 0);
        }

        /* ------------------------------------------------------------------------ */

        void
        XTKParameterList(Module_Parameter_Lists &aParameterLists)
        {
            aParameterLists.set("decompose", true);
            aParameterLists.set("decomposition_type", "conformal");
            aParameterLists.set("enrich_mesh_indices", "0");
            aParameterLists.set("ghost_stab", true);
            aParameterLists.set("multigrid", false);
            aParameterLists.set("verbose", false);
            aParameterLists.set("high_to_low_dbl_side_sets", true);
            aParameterLists.set("print_enriched_ig_mesh", false);
            aParameterLists.set("exodus_output_XTK_ig_mesh", true);
            // aParameterLists.set("only_generate_xtk_temp", true);
        }

        /* ------------------------------------------------------------------------ */

        void
        GENParameterList(Module_Parameter_Lists &aParameterLists)
        {

            aParameterLists.set("IQI_types", "IQIPerimeter");

            aParameterLists.set("number_of_phases", 4);
            aParameterLists.set("phase_function_name", "get_phase_index");
            // aParameterLists.set("output_mesh_file", tGENOutputFile);
            aParameterLists.set("time_offset", 10.0);

            // Height design variable
            Design_Variable tHeightDV(-0.1, 0.0, 0.1);
            Design_Variable tThetaDV(0.0, M_PI / 6.0, M_PI / 2.0); // High Error
            // Design_Variable tThetaDV(0.0, 0.273181969877373, M_PI / 2.0); // High Error
            // Design_Variable tThetaDV(0.0, 0.462307949023247, M_PI / 2.0); // Low Error

            // Fixed bottom line
            aParameterLists(GEN::GEOMETRIES).add_parameter_list(gen::Field_Type::USER_DEFINED);
            aParameterLists.set("name", "Bottom_Line");
            aParameterLists.set("field_function_name", "Bottom_Line");
            aParameterLists.set("sensitivity_function_name", "dLine_dHeight");
            aParameterLists.insert("theta", tThetaDV);
            aParameterLists.set("compute_sensitivities_along_edges", tUseConsistentSens);
            aParameterLists.set("use_multilinear_interpolation", tUseMultilinearInterpolation);

            // Fixed left line
            aParameterLists(GEN::GEOMETRIES).add_parameter_list(gen::Field_Type::USER_DEFINED);
            aParameterLists.set("name", "Left_Line");
            aParameterLists.set("field_function_name", "Left_Line");
            aParameterLists.set("sensitivity_function_name", "dLine_dHeight");
            aParameterLists.insert("theta", tThetaDV);
            aParameterLists.set("compute_sensitivities_along_edges", tUseConsistentSens);
            aParameterLists.set("use_multilinear_interpolation", tUseMultilinearInterpolation);

            // Fixed right line
            aParameterLists(GEN::GEOMETRIES).add_parameter_list(gen::Field_Type::USER_DEFINED);
            aParameterLists.set("name", "Right_Line");
            aParameterLists.set("field_function_name", "Right_Line");
            aParameterLists.set("sensitivity_function_name", "dLine_dHeight");
            aParameterLists.insert("theta", tThetaDV);
            aParameterLists.set("compute_sensitivities_along_edges", tUseConsistentSens);
            aParameterLists.set("use_multilinear_interpolation", tUseMultilinearInterpolation);

            // Movable top line - important this is LAST!
            aParameterLists(GEN::GEOMETRIES).add_parameter_list(gen::Field_Type::USER_DEFINED);
            aParameterLists.set("name", "Top_Line");
            aParameterLists.set("field_function_name", "Top_Line");
            aParameterLists.set("sensitivity_function_name", "dTop_Line_dHeight");
            aParameterLists.set("spatial_gradient_function_name", "dTop_Line_dX");
            aParameterLists.insert("theta", tThetaDV);
            aParameterLists.insert("height", tHeightDV);
            aParameterLists.set("compute_sensitivities_along_edges", tUseConsistentSens);
            aParameterLists.set("use_multilinear_interpolation", tUseMultilinearInterpolation);
        }

        /* ------------------------------------------------------------------------ */

        void
        FEMParameterList(Module_Parameter_Lists &aParameterLists)
        {

            //------------------------------------------------------------------------------
            if (tUsePhases)
            {
                aParameterLists(FEM::PHASES).add_parameter_list();
                aParameterLists.set("phase_name", "PhaseSolid");
                aParameterLists.set("phase_indices", "0");

                aParameterLists(FEM::PHASES).add_parameter_list();
                aParameterLists.set("phase_name", "PhaseVoid");
                aParameterLists.set("phase_indices", "1");

                aParameterLists(FEM::PHASES).add_parameter_list();
                aParameterLists.set("phase_name", "PhaseBottom");
                aParameterLists.set("phase_indices", "2");
            }
            else
            {
                aParameterLists.hack_for_legacy_fem();
            }

            //------------------------------------------------------------------------------
            // fill the property part of the parameter list

            // soild properties ------------------------------------------------------------

            aParameterLists(FEM::PROPERTIES).add_parameter_list();
            aParameterLists.set("property_name", "PropDensity");
            aParameterLists.set("function_parameters", tDensity);

            aParameterLists(FEM::PROPERTIES).add_parameter_list();
            aParameterLists.set("property_name", "PropCapacity");
            aParameterLists.set("function_parameters", tCapacity);

            // create  conductivity property
            aParameterLists(FEM::PROPERTIES).add_parameter_list();
            aParameterLists.set("property_name", "PropConductivity");
            aParameterLists.set("function_parameters", tConductivity);

            aParameterLists(FEM::PROPERTIES).add_parameter_list();
            aParameterLists.set("property_name", "PropInletTemp");
            aParameterLists.set("function_parameters", "0.0");

            aParameterLists(FEM::PROPERTIES).add_parameter_list();
            aParameterLists.set("property_name", "PropOutletFlux");
            aParameterLists.set("function_parameters", "1.0");

            // BC properties ---------------------------------------------------------------

            // create heat load property
            aParameterLists(FEM::PROPERTIES).add_parameter_list();
            aParameterLists.set("property_name", "PropVolumetricHeatFlux");
            aParameterLists.set("value_function", "Func_HeatLoad");

            //------------------------------------------------------------------------------
            // fill the constitutive model part of the parameter list

            //  CM --------------------------------------------------------------------
            // create  diffusion CM
            aParameterLists(FEM::CONSTITUTIVE_MODELS).add_parameter_list();
            aParameterLists.set("constitutive_name", "CMDiffusion");
            aParameterLists.set("constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO);
            aParameterLists.set("dof_dependencies", std::pair<std::string, std::string>("TEMP", "Temperature"));
            aParameterLists.set("properties",
                                "PropConductivity,Conductivity;"
                                "PropDensity     ,Density;"
                                "PropCapacity    ,HeatCapacity");
            if (tUsePhases)
            {
                aParameterLists.set("phase_name", "PhaseSolid");
            }

            //------------------------------------------------------------------------------
            // fill the stabilization parameter part of the parameter list

            // create Nitsche for  temperature
            aParameterLists(FEM::STABILIZATION).add_parameter_list();
            aParameterLists.set("stabilization_name", "SPNitscheT");
            aParameterLists.set("stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE);
            aParameterLists.set("function_parameters", tNitschePenThermal);
            aParameterLists.set("leader_properties", "PropConductivity,Material");
            if (tUsePhases)
            {
                aParameterLists.set("phase_name", "PhaseSolid");
            }

            // create ghost penalty  temperature
            aParameterLists(FEM::STABILIZATION).add_parameter_list();
            aParameterLists.set("stabilization_name", "SPGPTemp");
            aParameterLists.set("stabilization_type", fem::Stabilization_Type::GHOST_DISPL);
            aParameterLists.set("function_parameters", "0.05");
            aParameterLists.set("leader_properties", "PropConductivity,Material");
            if (tUsePhases)
            {
                aParameterLists.set("leader_phase_name", "PhaseSolid");
                aParameterLists.set("follower_phase_name", "PhaseSolid");
            }

            //------------------------------------------------------------------------------
            // fill the IWG part of the parameter list

            //  bulk IWGs -------------------------------------------------------------
            // diffusion
            aParameterLists(FEM::IWG).add_parameter_list();
            aParameterLists.set("IWG_name", "IWGDiffusionBulk");
            aParameterLists.set("IWG_type", fem::IWG_Type::SPATIALDIFF_BULK);
            aParameterLists.set("dof_residual", "TEMP");
            aParameterLists.set("leader_properties", "PropVolumetricHeatFlux,Load");
            aParameterLists.set("leader_constitutive_models", "CMDiffusion,Diffusion");
            if (tUsePhases)
            {
                aParameterLists.set("phase_name", "PhaseSolid");
            }
            else
            {
                aParameterLists.set("mesh_set_names", tSolidPhase);
            }

            aParameterLists(FEM::IWG).add_parameter_list();
            aParameterLists.set("IWG_name", "IWGInletTemp");
            aParameterLists.set("IWG_bulk_type", fem::Element_Type::SIDESET);
            aParameterLists.set("IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE);
            aParameterLists.set("dof_residual", "TEMP");
            aParameterLists.set("leader_properties", "PropInletTemp,Dirichlet");
            aParameterLists.set("leader_constitutive_models", "CMDiffusion,Diffusion");
            aParameterLists.set("stabilization_parameters", "SPNitscheT,DirichletNitsche");
            if (tUsePhases)
            {
                aParameterLists.set("leader_phase_name", "PhaseSolid");
                aParameterLists.set("neighbor_phases", "PhaseBottom");
            }
            else
            {
                aParameterLists.set("mesh_set_names", "iside_b0_0_b1_2");
            }

            aParameterLists(FEM::IWG).add_parameter_list();
            aParameterLists.set("IWG_name", "IWGOutletFlux");
            // aParameterLists.set("IWG_bulk_type", fem::Element_Type::SIDESET);
            aParameterLists.set("IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN);
            aParameterLists.set("dof_residual", "TEMP");
            aParameterLists.set("leader_dof_dependencies", "TEMP");
            aParameterLists.set("leader_properties", "PropOutletFlux,Neumann");
            if (tUsePhases)
            {
                aParameterLists.set("leader_phase_name", "PhaseSolid");
                aParameterLists.set("neighbor_phases", "PhaseVoid");
            }
            else
            {
                aParameterLists.set("mesh_set_names", tInterface);
            }

            if (tUseGhost)
            {
                // ghost  temperature
                aParameterLists(FEM::IWG).add_parameter_list();
                aParameterLists.set("IWG_name", "IWGGPTemp");
                aParameterLists.set("IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD);
                aParameterLists.set("dof_residual", "TEMP");
                aParameterLists.set("leader_dof_dependencies", "TEMP");
                aParameterLists.set("follower_dof_dependencies", "TEMP");
                aParameterLists.set("stabilization_parameters", "SPGPTemp,GhostSP");
                aParameterLists.set("ghost_order", 1);
                if (tUsePhases)
                {
                    aParameterLists.set("leader_phase_name", "PhaseSolid");
                    aParameterLists.set("follower_phase_name", "PhaseSolid");
                }
                else
                {
                    aParameterLists.set("mesh_set_names", tSolidPhase);
                }
            }

            //------------------------------------------------------------------------------
            // fill the IQI part of the parameter list

            aParameterLists(FEM::IQI).add_parameter_list();
            aParameterLists.set("IQI_name", "IQINodeX");
            aParameterLists.set("IQI_type", fem::IQI_Type::NODE_SENSITIVITY);
            aParameterLists.set("IQI_bulk_type", fem::Element_Type::SIDESET);
            aParameterLists.set("leader_dv_dependencies", "DXDP");
            aParameterLists.set("vectorial_field_index", 0);
            aParameterLists.set("function_parameters", "2");
            if (tUsePhases)
            {
                aParameterLists.set("phase_name", "PhaseSolid");
                aParameterLists.set("neighbor_phases", "PhaseVoid");
                aParameterLists.set("follower_phase_name", "PhaseVoid");
            }
            else
            {
                aParameterLists.set("mesh_set_names", tInterface);
            }

            aParameterLists(FEM::IQI).add_parameter_list();
            aParameterLists.set("IQI_name", "IQINodeY");
            aParameterLists.set("IQI_type", fem::IQI_Type::NODE_SENSITIVITY);
            aParameterLists.set("IQI_bulk_type", fem::Element_Type::SIDESET);
            aParameterLists.set("leader_dv_dependencies", "DXDP");
            aParameterLists.set("vectorial_field_index", 1);
            aParameterLists.set("function_parameters", "2");
            if (tUsePhases)
            {
                aParameterLists.set("phase_name", "PhaseSolid");
                aParameterLists.set("neighbor_phases", "PhaseVoid");
                aParameterLists.set("follower_phase_name", "PhaseVoid");
            }
            else
            {
                aParameterLists.set("mesh_set_names", tInterface);
            }

            // input thermal energy
            aParameterLists(FEM::IQI).add_parameter_list();
            aParameterLists.set("IQI_name", "IQIPerimeter");
            aParameterLists.set("IQI_bulk_type", fem::Element_Type::SIDESET);
            aParameterLists.set("IQI_type", (fem::IQI_Type::VOLUME));
            if (tUsePhases)
            {
                aParameterLists.set("phase_name", "PhaseSolid");
                aParameterLists.set("neighbor_phases", "PhaseVoid");
                aParameterLists.set("follower_phase_name", "PhaseVoid");
            }
            else
            {
                aParameterLists.set("mesh_set_names", tInterface);
            }

            //------------------------------------------------------------------------------
            // fill the computation part of the parameter list
            aParameterLists(FEM::COMPUTATION);
            aParameterLists.set("print_physics_model", false);

            aParameterLists.set("finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL);
            aParameterLists.set("finite_difference_perturbation_size", tFEMFdEpsilon);
            aParameterLists.set("finite_difference_perturbation_strategy", fem::Perturbation_Type::ABSOLUTE);
        }

        void
        SOLParameterList(Module_Parameter_Lists &aParameterLists)
        {

            aParameterLists(SOL::LINEAR_ALGORITHMS).add_parameter_list(sol::SolverType::AMESOS_IMPL);
#ifdef MORIS_USE_MUMPS
            aParameterLists.set("Solver_Type", "Amesos_Mumps");
#else
        aParameterLists.set("Solver_Type", "Amesos_Superludist");
#endif

            //------------------------------------------------------------------------------

            aParameterLists(SOL::LINEAR_SOLVERS).add_parameter_list();
            aParameterLists.set("DLA_Linear_solver_algorithms", "0");

            //------------------------------------------------------------------------------

            aParameterLists(SOL::NONLINEAR_ALGORITHMS).add_parameter_list(); // nonlinear algorithm index 0
            aParameterLists.set("NLA_Solver_Implementation", NLA::NonlinearSolverType::NEWTON_SOLVER);
            aParameterLists.set("NLA_Linear_solver", 0);
            aParameterLists.set("NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop);
            aParameterLists.set("NLA_relaxation_parameter", tNLA_relaxation_parameter);
            aParameterLists.set("NLA_max_iter", tNLA_max_iter);

            //------------------------------------------------------------------------------

            aParameterLists(SOL::NONLINEAR_SOLVERS).add_parameter_list(); // 1: thermal subproblem
            aParameterLists.set("NLA_Nonlinear_solver_algorithms", "0");  // set nonlinear algorithm with index 0
            aParameterLists.set("NLA_Solver_Implementation", NLA::NonlinearSolverType::NEWTON_SOLVER);
            aParameterLists.set("NLA_DofTypes", "TEMP");

            // ----------------------------------------------------------

            aParameterLists(SOL::TIME_SOLVER_ALGORITHMS).add_parameter_list();
            aParameterLists.set("TSA_Nonlinear_Solver", 0);             // using NLBGS for forward problem
            aParameterLists.set("TSA_Nonlinear_Sensitivity_Solver", 0); // using monlithic for sensitivity problem

            //------------------------------------------------------------------------------

            aParameterLists(SOL::TIME_SOLVERS).add_parameter_list();
            aParameterLists.set("TSA_DofTypes", "TEMP");
            aParameterLists.set("TSA_Initialize_Sol_Vec", "TEMP,0.0");
            aParameterLists.set("TSA_Output_Indices", "0");
            aParameterLists.set("TSA_Output_Criteria", "Output_Criterion");
            aParameterLists.set("TSA_time_level_per_type", "TEMP,1");

            //------------------------------------------------------------------------------

            aParameterLists(SOL::SOLVER_WAREHOUSE).set("SOL_save_operator_to_matlab", "Jacobian");
            aParameterLists.set("Sensitivity_Analysis_Type", sol::SensitivityAnalysisType::DIRECT);

            aParameterLists(SOL::PRECONDITIONERS).add_parameter_list(sol::PreconditionerType::NONE);
        }

        void
        MSIParameterList(Module_Parameter_Lists &aParameterLists)
        {
            aParameterLists.set("TEMP", 0);
        }

        void
        VISParameterList(Module_Parameter_Lists &aParameterLists)
        {
            aParameterLists.set("File_Name", std::pair<std::string, std::string>("./", tExoFile));
            aParameterLists.set("Mesh_Type", vis::VIS_Mesh_Type::STANDARD);
            aParameterLists.set("Set_Names", tAllPhases);
            aParameterLists.set("Field_Names", "PERIMETER,DX,DY");
            aParameterLists.set("Field_Type", "GLOBAL,NODAL,NODAL");
            aParameterLists.set("IQI_Names", "IQIPerimeter,IQINodeX,IQINodeY");
            aParameterLists.set("Analysis_Type", "SENSITIVITY,SENSITIVITY,SENSITIVITY");
            aParameterLists.set("Save_Frequency", 1);
            aParameterLists.set("Time_Offset", 10.0);
        }

        void
        MORISGENERALParameterList(Module_Parameter_Lists &aParameterLists)
        {
        }

        //------------------------------------------------------------------------------
    } // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
