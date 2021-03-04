/* SCHNAKENBERG REACTION_DIFFUSION SYSTEM with periodic boundary conditions
 *
 * Author: Aaditya Lakshmanan, University of Michigan, 2020
 */


// Include files
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/generic_linear_algebra.h>

// For reading parameter files 
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/parameter_handler.h>
#include <list>


#include <fstream>
#include <iostream>

// To generate a random integer
#include <experimental/random>


// Solvers
namespace LA
{
	#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
	!(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
	using namespace dealii::LinearAlgebraPETSc;
	#  define USE_PETSC_LA
	#elif defined(DEAL_II_WITH_TRILINOS)
	using namespace dealii::LinearAlgebraTrilinos;
	#else
	#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
	#endif
} // namespace LA

// Headers
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/distributed/tria.h>



#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>



#include <deal.II/numerics/data_out.h>



#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/utilities.h>


#include <deal.II/base/conditional_ostream.h>

#include <deal.II/base/index_set.h>



#include <deal.II/lac/sparsity_tools.h>



	
// The last step is as in all previous programs:
namespace RD_Schnakenberg
{
	using namespace dealii;
  
    // To read in parameters
    ParameterHandler     prm;
   
    // Parameter entries
    void declare_parameters (){
		prm.declare_entry ("gamma", "220.0",
                           Patterns::Double(),
                           "Gamma in the Schnakenberg model");
        prm.declare_entry ("apar", "0.2",
                           Patterns::Double(),
						   "a in the Schnakenberg model");
        prm.declare_entry ("bpar", "1.3",
							Patterns::Double(),
						   "b in the Schnakenberg model");
        prm.declare_entry ("diffusivity", "119.0",
                           Patterns::Double(),
                           "d in the Schnakenberg model");
	    prm.declare_entry ("initial perturbation", "0.005",
                            Patterns::Double(),
                           "Perturbation in the Schnakenberg model");
	    prm.declare_entry ("xlow", "0.0",
                           Patterns::Double(),
                           "x lower coordinate");
        prm.declare_entry ("xhigh", "1.0",
                           Patterns::Double(),
                           "x upper coordinate");	
        prm.declare_entry ("ylow", "0.0",
		                    Patterns::Double(),         
			               "y lower coordinate");
        prm.declare_entry ("yhigh", "0.01",
		                    Patterns::Double(),
			               "y upper coordinate");
        prm.declare_entry ("x subdivisions", "25",
							Patterns::Integer(),
							"x direction refinement");
        prm.declare_entry ("y subdivisions", "1",
                            Patterns::Integer(),
                          "y direction refinement");
	    prm.declare_entry ("time step", "0.00083333",
							Patterns::Double(),
							"time step");
        prm.declare_entry ("total time", "0.08333333",
						    Patterns::Double(),
                            "total time");
        prm.declare_entry ("skip frequency", "5",
							Patterns::Integer(),
							"skip frequency");				   
   
    }
  
	// Parameters for the Schnakenberg system
	double gamma, a_param, b_param, d_param ;
	double prtrb ;
	// Mesh parameters
    double min_x_coord, max_x_coord, min_y_coord, max_y_coord  ;
    unsigned int x_subdivs, y_subdivs ;
	double       time_step;
    double       total_time ;
	unsigned int output_skip;
	
	// Read from parameters file
	void parse_command_line (const int     argc,
							char *const *argv){
		if (argc < 2)
		{
     		std::cout << "Not enough input arguements" << std::endl ;

			exit (1);
		}
	    std::list<std::string> args;
			
		for (int i=0; i<argc; ++i){
			args.emplace_back(argv[i]);
		}	
		std::string parameter_file ;
		if (args.front() == std::string("mpirun")){
		    args.pop_front();
			args.pop_front();
			args.pop_front();
			args.pop_front();
			parameter_file = args.front ();
			args.pop_front();
			}
		else{ 
			args.pop_front();
			parameter_file = args.front ();
			args.pop_front();
		}
			
		prm.parse_input (parameter_file);			
			
		gamma = prm.get_double("gamma") ;
		a_param = prm.get_double("apar") ;
		b_param = prm.get_double("bpar") ;
		d_param = prm.get_double("diffusivity") ;
		prtrb = prm.get_double("initial perturbation") ;
		min_x_coord = prm.get_double("xlow") ;
		max_x_coord = prm.get_double("xhigh") ;
		min_y_coord = prm.get_double("ylow") ;
		max_y_coord = prm.get_double("yhigh") ;
		x_subdivs = prm.get_integer("x subdivisions") ;
		y_subdivs = prm.get_integer("y subdivisions") ;
		time_step = prm.get_double("time step") ;
		total_time = prm.get_double("total time") ;
		output_skip = prm.get_integer("skip frequency") ;
		
	}						
	  
	// Declaration of the problem class
	template <int dim>
	class Schnakenberg
	{
		public:
			Schnakenberg();
			void run();
	

		private:
			void create_mesh();
			void setup_system();
			void assemble_u();
			void assemble_v();
			void solve_u();
			void solve_v();
			void output_results(const unsigned int timestep_number) const;
	
		MPI_Comm mpi_communicator;

		parallel::distributed::Triangulation<dim> triangulation;
	
		FE_Q<dim>          fe;
		DoFHandler<dim>    dof_handler;
    
		IndexSet locally_owned_dofs;
		IndexSet locally_relevant_dofs; 
 	
		AffineConstraints<double> constraints;
    
	
		LA::MPI::SparseMatrix matrix_u;
		LA::MPI::SparseMatrix matrix_v;
		LA::MPI::Vector  locally_relevant_solution_u, locally_relevant_solution_v
						, locally_relevant_old_solution_u, locally_relevant_old_solution_v
						, system_rhs_u, system_rhs_v;
	
		ConditionalOStream pcout;
		TimerOutput        computing_timer;
	
		// Time information
		double   time ;
	};



// 	Equation data for initial values of U
	template <int dim>
	class InitialValuesU : public Function<dim>
	{
		public:
    
			virtual double value(const Point<dim> & p,
								const unsigned int component = 0) const override{    
	   
				unsigned int randnumint ;
				std::experimental::reseed () ;
				randnumint = std::experimental::randint(1,1000000) ;
				srand (randnumint) ;

				if(p[0] == min_x_coord || p[0] == max_x_coord || p[1] == min_y_coord || p[1] == max_y_coord )
					return (a_param + b_param) ;
				else 	
					return ((a_param + b_param)*( 1.0 + prtrb*static_cast <double> (rand()) / static_cast <double> (RAND_MAX) ));
		
//				return (a_param + b_param)*( 1. + ( prtrb * std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]) ) ) ;
			}
	};


//	Equation data for initial values of U
	template <int dim>
	class InitialValuesV : public Function<dim>
	{
		public:
			virtual double value(const Point<dim> & p,
								const unsigned int component = 0) const override{      
				unsigned int randnumint ;
				std::experimental::reseed () ;
				randnumint = std::experimental::randint(1,1000000) ;
				srand (randnumint) ; 
				if(p[0] == min_y_coord || p[0] == max_y_coord || p[1] == min_y_coord || p[1] == max_y_coord )
					return (b_param/(a_param+b_param)/(a_param+b_param)) ;
				else 
					return ((b_param/(a_param+b_param)/(a_param+b_param))*(1.0  + prtrb*static_cast <double> (rand()) / static_cast <double> (RAND_MAX) ));
	   
//	  			return b_param/(a_param + b_param)/(a_param + b_param)*( 1. + ( prtrb * std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]) ) ) ;
	   
			}
	};



// 	Problem class constructor to initialize member variables
	template <int dim>
	Schnakenberg<dim>::Schnakenberg()
    : mpi_communicator(MPI_COMM_WORLD)
	, triangulation(mpi_communicator)
	, fe(1)
    , dof_handler(triangulation)
	, pcout(std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
	, computing_timer(mpi_communicator,
					  pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)	
	, time(time_step)
   {}


//   Mesh and information about periodicity constraints 
	template <int dim>
	void Schnakenberg<dim>::create_mesh(){
		
	    TimerOutput::Scope t(computing_timer, "setup");
  
//  	GridGenerator::hyper_cube(triangulation, min_coord, max_coord, true);
		std::vector<unsigned int> subdivs_list ;
		subdivs_list.push_back(x_subdivs) ;
		subdivs_list.push_back(y_subdivs) ;
		GridGenerator::subdivided_hyper_rectangle(triangulation, subdivs_list,
												  Point<2>(min_x_coord,min_y_coord),
												  Point<2>(max_x_coord,max_y_coord),true);
	  
		std::vector<GridTools::PeriodicFacePair<
					typename parallel::distributed::Triangulation<dim>::cell_iterator>>
					periodicity_vector;
	  
	  /*for (auto &face : triangulation.active_face_iterators())
		if (face->at_boundary())
			if (face->center()[0] == min_coord)
				face->set_boundary_id (1);
	  
      for (auto &face : triangulation.active_face_iterators())
		if (face->at_boundary())
			if (face->center()[0] == max_coord)
				face->set_boundary_id (2);
	  
      for (auto &face : triangulation.active_face_iterators())
		if (face->at_boundary())
			if (face->center()[1] == min_coord)
				face->set_boundary_id (3);
				
			

      for (auto &face : triangulation.active_face_iterators())
		if (face->at_boundary())
			if (face->center()[1] == max_coord)
				face->set_boundary_id (4); */
							
	  		
		GridTools::collect_periodic_faces(triangulation,
										  0,
									      1,
										  0,
										  periodicity_vector);
										
	  									
	  
		GridTools::collect_periodic_faces(triangulation,
										  2,
										  3,
										  1,
										  periodicity_vector); 	
	
	  
	  
		triangulation.add_periodicity(periodicity_vector);
			  
	}

  
//  Setting up the problem - constraints, initializing solution variables
	template <int dim>
	void Schnakenberg<dim>::setup_system(){ 
    
	
		dof_handler.distribute_dofs(fe);	
	
		locally_owned_dofs = dof_handler.locally_owned_dofs();

		DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
      
	
		constraints.clear();
		constraints.reinit(locally_relevant_dofs);
    
		std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
					periodicity_vector;
		
		GridTools::collect_periodic_faces(dof_handler,
										  0,
										  1,
										  0,
										  periodicity_vector);    
					
					
		GridTools::collect_periodic_faces(dof_handler,
										  2,
										  3,
										  1,
										  periodicity_vector);    
									
		DoFTools::make_periodicity_constraints<DoFHandler<dim>>(periodicity_vector,
																constraints);				
    
		constraints.close();
    	
		DynamicSparsityPattern dsp(locally_relevant_dofs); 
		DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
		SparsityTools::distribute_sparsity_pattern(dsp,
												   dof_handler.locally_owned_dofs(),
                                                   mpi_communicator,
                                                   locally_relevant_dofs); 
    
		locally_relevant_solution_u.reinit(locally_owned_dofs,
										   locally_relevant_dofs,
                                           mpi_communicator);
		locally_relevant_old_solution_u.reinit(locally_owned_dofs,
	                                           locally_relevant_dofs,
                                               mpi_communicator);
		locally_relevant_solution_v.reinit(locally_owned_dofs,
	                                   locally_relevant_dofs, 
                                       mpi_communicator);
		locally_relevant_old_solution_v.reinit(locally_owned_dofs,
										   locally_relevant_dofs,
                                           mpi_communicator);
										   
		system_rhs_u.reinit(locally_owned_dofs, mpi_communicator);
		system_rhs_v.reinit(locally_owned_dofs, mpi_communicator);
	
		matrix_u.reinit(locally_owned_dofs,
                        locally_owned_dofs,
                        dsp,
                        mpi_communicator);
	
		matrix_v.reinit(locally_owned_dofs,
                        locally_owned_dofs,
                        dsp,
                        mpi_communicator);  	
    
	
	}

// 	Assembly for the U equation
	template <int dim>
	void Schnakenberg<dim>::assemble_u(){
   	  
		TimerOutput::Scope t(computing_timer, "assembly");
  
		const QGauss<dim> quadrature_formula(fe.degree + 1);
		FEValues<dim> fe_values(fe,
								quadrature_formula,
								update_values | update_gradients |
								update_quadrature_points | update_JxW_values);
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		const unsigned int n_q_points    = quadrature_formula.size();
  
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		Vector<double>     cell_rhs(dofs_per_cell);
  

		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		std::vector<double>                  solution_u_gp(n_q_points);
		std::vector<double>                  solution_v_gp(n_q_points);
  
   
		for (const auto &cell : dof_handler.active_cell_iterators()){
			if (cell->is_locally_owned())
			{
				cell_matrix = 0.;
				cell_rhs    = 0.;
				fe_values.reinit(cell);
				fe_values.get_function_values(locally_relevant_old_solution_u, solution_u_gp);
				fe_values.get_function_values(locally_relevant_old_solution_v, solution_v_gp);
		
				for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
				{
					// U equation matrix 
					for (unsigned int i = 0; i < dofs_per_cell; ++i)
					{
						for (unsigned int j = 0; j < dofs_per_cell; ++j){
							cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
												fe_values.shape_grad(j, q_point) *
												fe_values.JxW(q_point)*time_step;
									   
							/*cell_matrix(i, j) += fe_values.shape_value(i, q_point) *
											      fe_values.shape_value(j, q_point) *
												  fe_values.JxW(q_point)*(1.0 + gamma*time_step);	*/
							cell_matrix(i, j) += fe_values.shape_value(i, q_point) *
												fe_values.shape_value(j, q_point) *
												fe_values.JxW(q_point) ;					   
						}
				
						cell_rhs(i) += gamma*solution_u_gp[q_point]*solution_u_gp[q_point]*solution_v_gp[q_point]
										*fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)*time_step;
						cell_rhs(i) += solution_u_gp[q_point]
										*fe_values.shape_value(i, q_point) * fe_values.JxW(q_point);
						cell_rhs(i) += -1.0*solution_u_gp[q_point]
										*fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)*gamma*time_step;				
				
						cell_rhs(i) += gamma*a_param*fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)*time_step;
					}
			  
				}
		  
		  
				cell->get_dof_indices(local_dof_indices);
		
				constraints.distribute_local_to_global(cell_matrix,
														cell_rhs,
														local_dof_indices,
														matrix_u,
														system_rhs_u);
			}
		}	
  
		matrix_u.compress(VectorOperation::add);
		system_rhs_u.compress(VectorOperation::add);
	}


	// Solve for U field 
	template <int dim>
	void Schnakenberg<dim>::solve_u(){
		
		TimerOutput::Scope t(computing_timer, "solve");
		LA::MPI::Vector    completely_distributed_solution(locally_owned_dofs,
														mpi_communicator);
		SolverControl solver_control(dof_handler.n_dofs(), 1e-12);
		#ifdef USE_PETSC_LA
		LA::SolverCG solver(solver_control, mpi_communicator);
		#else
		LA::SolverCG solver(solver_control);
		#endif
		LA::MPI::PreconditionAMG preconditioner;
		LA::MPI::PreconditionAMG::AdditionalData data;
		#ifdef USE_PETSC_LA
		data.symmetric_operator = true;
		#else
		/* Trilinos defaults are good */
		#endif
		preconditioner.initialize(matrix_u, data);
		solver.solve(matrix_u,
					completely_distributed_solution,
					system_rhs_u,
					preconditioner);
	
    
	
		pcout << "   Solved in " << solver_control.last_step() << " iterations."
			  << std::endl;
		constraints.distribute(completely_distributed_solution);
		locally_relevant_solution_u = completely_distributed_solution;
		locally_relevant_old_solution_u = completely_distributed_solution;
	
		matrix_u = 0. ;	
		system_rhs_u = 0. ;
	}

	// Assembly for the V equation
	template <int dim>
	void Schnakenberg<dim>::assemble_v(){
	  
		TimerOutput::Scope t(computing_timer, "assembly");
		const QGauss<dim> quadrature_formula(fe.degree + 1);
		FEValues<dim> fe_values(fe,
								quadrature_formula,
								update_values | update_gradients |
								update_quadrature_points | update_JxW_values);
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		const unsigned int n_q_points    = quadrature_formula.size();
  
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		Vector<double>     cell_rhs(dofs_per_cell);
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		std::vector<double>                  solution_u_gp(n_q_points);
		std::vector<double>                  solution_v_gp(n_q_points);
  
  
		for (const auto &cell : dof_handler.active_cell_iterators()){
			if (cell->is_locally_owned())
			{
				cell_matrix = 0.;
				cell_rhs    = 0.;
				fe_values.reinit(cell);
				fe_values.get_function_values(locally_relevant_old_solution_u, solution_u_gp);
				fe_values.get_function_values(locally_relevant_old_solution_v, solution_v_gp);
		
				for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
				{
					// V equation matrix 
					for (unsigned int i = 0; i < dofs_per_cell; ++i)
					{
						for (unsigned int j = 0; j < dofs_per_cell; ++j){
							cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
												fe_values.shape_grad(j, q_point) *
												fe_values.JxW(q_point)*time_step*d_param;
									   
							cell_matrix(i, j) += fe_values.shape_value(i, q_point) *
												fe_values.shape_value(j, q_point) *
												fe_values.JxW(q_point);	
						}									   
						// V equation RHS			   
						cell_rhs(i) += -1.0*gamma*solution_u_gp[q_point]*solution_u_gp[q_point]*solution_v_gp[q_point]
										*fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)*time_step;
								
						cell_rhs(i) += solution_v_gp[q_point]
									*fe_values.shape_value(i, q_point) * fe_values.JxW(q_point);				
				
						cell_rhs(i) += gamma*b_param*fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)*time_step;
					}
			    }
		  
		  
				cell->get_dof_indices(local_dof_indices);
		
				constraints.distribute_local_to_global(cell_matrix,
														cell_rhs,
														local_dof_indices,
														matrix_v,
														system_rhs_v);
			}
		}
		matrix_v.compress(VectorOperation::add);
		system_rhs_v.compress(VectorOperation::add);
	}

	// Solve for the V field
	template <int dim>
	void Schnakenberg<dim>::solve_v()
	{
		TimerOutput::Scope t(computing_timer, "solve");
		LA::MPI::Vector    completely_distributed_solution(locally_owned_dofs,
														mpi_communicator);
		SolverControl solver_control(dof_handler.n_dofs(), 1e-12);
		#ifdef USE_PETSC_LA
		LA::SolverCG solver(solver_control, mpi_communicator);
		#else
		LA::SolverCG solver(solver_control);
		#endif
		LA::MPI::PreconditionAMG preconditioner;
		LA::MPI::PreconditionAMG::AdditionalData data;
		#ifdef USE_PETSC_LA
		data.symmetric_operator = true;
		#else
		/* Trilinos defaults are good */
		#endif
		preconditioner.initialize(matrix_v, data);
		solver.solve(matrix_v,
					completely_distributed_solution,
					system_rhs_v,
					preconditioner);
		pcout << "   Solved in " << solver_control.last_step() << " iterations."
			 << std::endl;
		constraints.distribute(completely_distributed_solution);
		locally_relevant_solution_v = completely_distributed_solution;
		locally_relevant_old_solution_v = completely_distributed_solution;
		matrix_v = 0. ;
		system_rhs_v = 0. ;
	}
    

	// Output results 
	template <int dim>
	void Schnakenberg<dim>::output_results(const unsigned int timestep_number) const
	{
		DataOut<dim> data_out;

		data_out.attach_dof_handler(dof_handler);
		data_out.add_data_vector(locally_relevant_solution_u, "U");
		data_out.add_data_vector(locally_relevant_solution_v, "V");

		Vector<float> subdomain(triangulation.n_active_cells());
		for (unsigned int i = 0; i < subdomain.size(); ++i){
			subdomain(i) = triangulation.locally_owned_subdomain();
		}		
		data_out.add_data_vector(subdomain, "subdomain");
		data_out.build_patches();
		data_out.write_vtu_with_pvtu_record(
				"./", "solution",timestep_number, mpi_communicator, 4, 8);
	}



	// Run method
	template <int dim>
	void Schnakenberg<dim>::run()
	{
	  
		pcout << "Running with "
		#ifdef USE_PETSC_LA
			  << "PETSc"
		#else
			  << "Trilinos"
		#endif
			  << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
              << " MPI rank(s)..." << std::endl;
	
		create_mesh() ;	
		setup_system();
		pcout << "Setup complete" ;
	
		LA::MPI::Vector    completely_distributed_solution_u(locally_owned_dofs,
														mpi_communicator);
		LA::MPI::Vector    completely_distributed_solution_v(locally_owned_dofs,
														mpi_communicator);												
		// Initialize U and V fields
		VectorTools::project(dof_handler,
							constraints,
							QGauss<dim>(fe.degree + 1),
							InitialValuesU<dim>(),
							completely_distributed_solution_u);
				 
		VectorTools::project(dof_handler,
							constraints,
							QGauss<dim>(fe.degree + 1),
							InitialValuesV<dim>(),
							completely_distributed_solution_v);
	
		locally_relevant_solution_u = completely_distributed_solution_u;
		locally_relevant_old_solution_u = completely_distributed_solution_u;
		locally_relevant_solution_v = completely_distributed_solution_v;
		locally_relevant_old_solution_v = completely_distributed_solution_v;
	
		pcout << "Applied initial condition" << std::endl ;
	
		   
		pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << 2*dof_handler.n_dofs()
              << std::endl;
	
    
	
		
		// Output initial value on mesh
		//   output_results(0);

    
		unsigned int timestep_number = 1;	
		//   Loop over time steps
		for (; time <= total_time; time += time_step, ++timestep_number)
		{
		
			pcout << "Time step " << timestep_number << " at t=" << time
				  << std::endl;
        
		
			// Assemble U and V equations
			assemble_u() ; 
			assemble_v() ;
		
			// Solve for U and V
			solve_u();
			solve_v();
		
			if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
			{
				TimerOutput::Scope t(computing_timer, "output");
            }

//        Output results based on time step skipping criterion
	  
			if (timestep_number % output_skip == 0)
		          	output_results(timestep_number);
		
		// Run time information
			computing_timer.print_summary();
			computing_timer.reset();
		
			pcout << std::endl;
    
		}
    // Only to write the output in the last timestep
		output_results(timestep_number-1) ;
	}
} // namespace end


// Main function
int main(int argc, char *argv[])
{
	try
    {
		using namespace dealii; 	 
		using namespace RD_Schnakenberg;
	  
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
	  
		declare_parameters ();
		parse_command_line (argc, argv);
	  
		Schnakenberg<2> schnakenberg_solver;
		schnakenberg_solver.run();
	}
	catch (std::exception &exc)
    {
		std::cerr << std::endl
					<< std::endl
					<< "----------------------------------------------------"
					<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
					<< exc.what() << std::endl
					<< "Aborting!" << std::endl
					<< "----------------------------------------------------"
					<< std::endl;

		return 1;
    }
	catch (...)
    {
		std::cerr << std::endl
					<< std::endl
					<< "----------------------------------------------------"
					<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
					<< "Aborting!" << std::endl
					<< "----------------------------------------------------"
					<< std::endl;
		return 1;
    }

	return 0;
}
