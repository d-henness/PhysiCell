/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

Cell_Definition* pImmuneCell;

void immune_cell_update_phenotype( Cell* pCell, Phenotype& phenotype, double dt ){
  static int start_phase_index = phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
  static int end_phase_index = phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
  static int apoptosis_index = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);


  // change amount of compound_x secreted
  if (pCell->phenotype.secretion.secretion_rates[0] - 1.0 < 0.0){
    pCell->phenotype.secretion.secretion_rates[0] = 0;
  }
  else{
    pCell->phenotype.secretion.secretion_rates[0] -= 1.0;
  }

  if(!pCell->phenotype.death.dead){
    // change reproduce rates
    double delta = (pCell->custom_data["replication_bonus"] / pCell->custom_data["base_cycle_length"]) - (1.0 / pCell->custom_data["base_cycle_length"]);
    if((phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) - (delta / 10.0)) < (1.0 / pCell->custom_data["base_cycle_length"])){
      phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = 1.0 / pCell->custom_data["base_cycle_length"];
    }
    else{
      phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) -= delta / 10.0;
    }
    // change death rates
  //  delta = (pCell->custom_data["max_death_rate"] - pCell->custom_data["base_apop_rate"]) / 10.0;
  //  if((pCell->phenotype.death.rates[apoptosis_index] + delta) < pCell->custom_data["max_death_rate"]){
  //    pCell->phenotype.death.rates[apoptosis_index] += delta;
  //  }

    // use a model based on a half life curve to update the rate of apoptosis
    pCell->phenotype.death.rates[apoptosis_index] = pCell->custom_data["base_apop_rate"] / std::exp(-pCell->custom_data["apop_rate_constant"] * pCell->custom_data["min_since_last_kill_attempt"]);
  }

  pCell->custom_data["min_since_last_kill_attempt"] += parameters.doubles("phenotype_update_time");
}

void create_immune_cell_type( void )
{
	pImmuneCell = find_cell_definition( "immune cell" );

	static int immuno_ID = microenvironment.find_density_index( "compound_x" );

	// reduce o2 uptake

	pImmuneCell->phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles("immune_relative_adhesion");
	pImmuneCell->phenotype.mechanics.cell_cell_repulsion_strength *= parameters.doubles("immune_relative_repulsion");

	// figure out mechanics parameters

	pImmuneCell->phenotype.mechanics.relative_maximum_attachment_distance = pImmuneCell->custom_data["max_attachment_distance"] / pImmuneCell->phenotype.geometry.radius;

	pImmuneCell->phenotype.mechanics.attachment_elastic_constant = pImmuneCell->custom_data["elastic_coefficient"];

	pImmuneCell->phenotype.mechanics.relative_detachment_distance = pImmuneCell->custom_data["max_attachment_distance" ] / pImmuneCell->phenotype.geometry.radius;

	// set functions

	pImmuneCell->functions.update_phenotype = immune_cell_update_phenotype;
	pImmuneCell->functions.custom_cell_rule = immune_cell_rule;
	pImmuneCell->functions.update_migration_bias = immune_cell_motility;
	pImmuneCell->functions.contact_function = adhesion_contact_function;

	// set custom data values

	return;
}

void create_cell_types( void )
{
	// set the random seed
	SeedRandom( parameters.ints("random_seed") );

	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	initialize_default_cell_definition();
	static int immuno_ID = microenvironment.find_density_index( "compound_x" ); // 1

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/
	cell_defaults.phenotype.mechanics.relative_maximum_attachment_distance =
		cell_defaults.custom_data["max_attachment_distance"] / cell_defaults.phenotype.geometry.radius;

	cell_defaults.phenotype.mechanics.relative_detachment_distance
		= cell_defaults.custom_data["max_attachment_distance"] / cell_defaults.phenotype.geometry.radius ;

	cell_defaults.phenotype.mechanics.attachment_elastic_constant
		= cell_defaults.custom_data[ "elastic_coefficient" ];

	//cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_and_immune_stimulation;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = adhesion_contact_function;
	cell_defaults.functions.update_migration_bias = NULL;

	cell_defaults.functions.update_phenotype = NULL;//phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	create_immune_cell_type();

	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	return;
}

void setup_microenvironment( void )
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

void introduce_immune_cells( void )
{
	double tumor_radius = -9e9; // 250.0;
	double temp_radius = 0.0;

	// for the loop, deal with the (faster) norm squared
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		temp_radius = norm_squared( (*all_cells)[i]->position );
		if( temp_radius > tumor_radius )
		{ tumor_radius = temp_radius; }
	}
	// now square root to get to radius
	tumor_radius = sqrt( tumor_radius );

	// if this goes wackadoodle, choose 250
	if( tumor_radius < 250.0 )
	{ tumor_radius = 250.0; }

	std::cout << "current tumor radius: " << tumor_radius << std::endl;

	// now seed immune cells

	int number_of_immune_cells =
		parameters.ints("number_of_immune_cells"); // 7500; // 100; // 40;
	double radius_inner = 0.0; //tumor_radius + parameters.doubles("initial_min_immune_distance_from_tumor");// 30.0; // 75 // 50;
	double radius_outer = 500.0;// radius_inner + parameters.doubles("thickness_of_immune_seeding_region"); // 75.0; // 100; // 1000 - 50.0;

	double mean_radius = 0.5*(radius_inner + radius_outer);
	double std_radius = 0.33*( radius_outer-radius_inner)/2.0;

	for( int i=0 ;i < number_of_immune_cells ; i++ )
	{
		double theta = UniformRandom() * 6.283185307179586476925286766559;
		double phi = acos( 2.0*UniformRandom() - 1.0 );

		double radius = NormalRandom( mean_radius, std_radius );

		Cell* pCell = create_cell( *pImmuneCell );
		pCell->assign_position( radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), 0.0);
	}

	return;
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	// create some of each type of cell

	Cell* pC;

  Cell_Definition* pCD = cell_definitions_by_index[0];
  std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
  for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
  {
    std::vector<double> position = {0,0,0};
    position[0] = (UniformRandom() * 100) - 50;
    position[1] = (UniformRandom() * 100) - 50;
    position[2] = 0.0;

    pC = create_cell( *pCD );
    pC->assign_position( position );
//    pC->phenotype.secretion.saturation_densities[0] = pC->custom_data["base_compound_x_secretion_target"] * UniformRandom();
//    pC->custom_data["current_compound_x_secretion_target"] = pC->phenotype.secretion.saturation_densities[0];
  }
	std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();

	return;
}


std::vector<std::string> my_coloring_function( Cell* pCell ){
	// immune are black
	std::vector< std::string > output( 4, "black" );

	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
    return output;
	}

	if( pCell->type == 1 )
	{
		output[0] = "lime";
		output[1] = "lime";
		output[2] = "green";
    return output;
	}

	// if I'm under attack, color me
	if( pCell->state.attached_cells.size() > 0 )
	{
		output[0] = "darkcyan"; // orangered // "purple"; // 128,0,128
		output[1] = "black"; // "magenta";
		output[2] = "cyan"; // "magenta"; //255,0,255
		return output;
	}

	// live cells are green, but shaded by oncoprotein value
	if( pCell->phenotype.death.dead == false )
	{
		int compound_x_saturation = (int) (pCell->phenotype.secretion.saturation_densities[0] / pCell->custom_data["base_compound_x_secretion_target"] * 255.0 );
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", compound_x_saturation, compound_x_saturation, 255-compound_x_saturation );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );

		return output;
	}

	// if not, dead colors


	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling ||
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed ||
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}

	return output;
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt ){
	if(phenotype.death.dead) return;

	static int start_phase_index;
	static int end_phase_index;
	static int apoptosis_index;
	static int compound_x_index = pCell->get_microenvironment()->find_density_index("compound_x");

  start_phase_index = phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
  end_phase_index = phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
  apoptosis_index = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);

	double pCompound_x = (pCell->nearest_density_vector())[compound_x_index]; // PhysiCell_constants::oxygen_index];
	int n = pCell->phenotype.cycle.current_phase_index();

	// this multiplier is for linear interpolation of the compound_x value
	double multiplier = 1.0;
	if(pCompound_x > pCell->custom_data["compound_x_thresh"]){
		multiplier = pCompound_x;
	}
	// now, update the appropriate cycle transition rate

	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier / pCell->custom_data["base_cycle_length"];

  // update cell based on pressure
  double pressure = pCell->state.simple_pressure;
  pCell->custom_data["pressure"] = pressure;
	// Update necrosis rate
	multiplier = 1.0;
	if(pressure > pCell->custom_data["pressure_thresh"])
	{
		multiplier = std::pow(pressure, 2);
	}

	// now, update the necrosis rate

	pCell->phenotype.death.rates[apoptosis_index] = multiplier * pCell->custom_data["base_apop_rate"];
  return;
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt ){
  return;
}

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt ){
  return;
}

void immune_cell_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if attached, biased motility towards director chemoattractant
	// otherwise, biased motility towards cargo chemoattractant

	static int immune_factor_index = microenvironment.find_density_index( "compound_x" );

	// if not docked, attempt biased chemotaxis
	if( pCell->state.attached_cells.size() == 0 )
	{
		phenotype.motility.is_motile = true;

		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(immune_factor_index);
		normalize( &( phenotype.motility.migration_bias_direction ) );
	}
	else
	{
		phenotype.motility.is_motile = false;
	}

	return;
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container();
	int i = 0;
	while( i < nearby.size() )
	{
		// don't try to kill yourself
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++;
	}

	return NULL;
}

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	static int compound_x_i = pTarget->custom_data.find_variable_index( "current_compound_x_secretion_target" );
	static int attach_rate_i = pAttacker->custom_data.find_variable_index( "attachment_rate" );

	double compound_x_secretion_saturation =
		pAttacker->custom_data["compound_x_secretion_saturation"];
	double compound_x_secretion_threshold =
		pAttacker->custom_data["compound_x_secretion_threshold"];
	double compound_x_secretion_difference = compound_x_secretion_saturation - compound_x_secretion_threshold;

	double max_attachment_distance =
		pAttacker->custom_data["max_attachment_distance"];
	double min_attachment_distance =
		pAttacker->custom_data["min_attachment_distance"];
	double attachment_difference = max_attachment_distance - min_attachment_distance;

  // do not attach to immune cells
  if (pTarget->type_name == "immune cell"){
    return false;
  }

	if( pTarget->custom_data[compound_x_i] > compound_x_secretion_threshold && pTarget->phenotype.death.dead == false )
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement );
		if( distance_scale > max_attachment_distance )
		{ return false; }

		double scale = pTarget->custom_data[compound_x_i];
		scale -= compound_x_secretion_threshold;
		scale /= compound_x_secretion_difference;
		if( scale > 1.0 )
		{ scale = 1.0; }

    // added bc all cells are the same
    scale = 1.0;
    // remove when done

		distance_scale *= -1.0;
		distance_scale += max_attachment_distance;
		distance_scale /= attachment_difference;
		if( distance_scale > 1.0 )
		{ distance_scale = 1.0; }

		if( UniformRandom() < pAttacker->custom_data[attach_rate_i] * scale * dt * distance_scale )
		{
//			std::cout << "\t attach!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl;
			attach_cells( pAttacker, pTarget );
		}

		return true;
	}

	return false;
}

bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt )
{
	static int oncoprotein_i = pTarget->custom_data.find_variable_index( "oncoprotein" );
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );
	static int kill_rate_index = pAttacker->custom_data.find_variable_index( "kill_rate" );

	double oncoprotein_saturation = pAttacker->custom_data["oncoprotein_saturation"]; // 2.0;
	double oncoprotein_threshold = pAttacker->custom_data["oncoprotein_threshold"]; // 0.5; // 0.1;
	double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;

	// new
	if( pTarget->custom_data[oncoprotein_i] < oncoprotein_threshold )
	{ return false; }

	// new
	double scale = pTarget->custom_data[oncoprotein_i];
	scale -= oncoprotein_threshold;
	scale /= oncoprotein_difference;
	if( scale > 1.0 )
	{ scale = 1.0; }

  // added bc all cells are the same
  scale = 1.0;
  // remove when done

	if( UniformRandom() < pAttacker->custom_data[kill_rate_index] * scale * dt )
	{
//		std::cout << "\t\t kill!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl;
//		update secretion_rates to simulate release of chemicals that attract other immune cells
    pAttacker->phenotype.secretion.secretion_rates[0] = 10.0;
    pAttacker->phenotype.secretion.saturation_densities[0] = 10.0;

    // briefly increase the repliation rate
    static int start_phase_index = pAttacker->phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
    static int end_phase_index = pAttacker->phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
    pAttacker->phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = pAttacker->custom_data["replication_bonus"] / pAttacker->custom_data["base_cycle_length"];

    // briefly return death_rate to normal
    static int apoptosis_index = pAttacker->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    pAttacker->phenotype.death.rates[apoptosis_index] = pAttacker->custom_data["base_apop_rate"];

		return true;
	}
	return false;
}

bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );

	// if the Target cell is already dead, don't bother!
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	return true;
}

void immune_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment_lifetime" );

	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions,
		// since those are part of mechanics.

		// Let's just fully disable now.
		pCell->functions.custom_cell_rule = NULL;
		return;
	}

	// if I'm docked
	if( pCell->state.number_of_attached_cells() > 0 )
	{
		// attempt to kill my attached cell

		bool detach_me = false;

    //update the kill timer
    pCell->custom_data["min_since_last_kill_attempt"] = 0.0;

		if( immune_cell_attempt_apoptosis( pCell, pCell->state.attached_cells[0], dt ) )
		{
			immune_cell_trigger_apoptosis( pCell, pCell->state.attached_cells[0] );
			detach_me = true;
		}

		// decide whether to detach

		if( UniformRandom() < dt / ( pCell->custom_data[attach_lifetime_i] + 1e-15 ) )
		{ detach_me = true; }

		// if I dettach, resume motile behavior

		if( detach_me )
		{
			detach_cells( pCell, pCell->state.attached_cells[0] );
			phenotype.motility.is_motile = true;
		}
		return;
	}

	// I'm not docked, look for cells nearby and try to docked

	// if this returns non-NULL, we're now attached to a cell
	if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
	{
		// set motility off
		phenotype.motility.is_motile = false;
		return;
	}
	phenotype.motility.is_motile = true;

	return;
}

void adhesion_contact_function( Cell* pActingOn, Phenotype& pao, Cell* pAttachedTo, Phenotype& pat , double dt )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position;

	static double max_elastic_displacement = pao.geometry.radius * pao.mechanics.relative_detachment_distance;
	static double max_displacement_squared = max_elastic_displacement*max_elastic_displacement;

	// detach cells if too far apart

	if( norm_squared( displacement ) > max_displacement_squared )
	{
		detach_cells( pActingOn , pAttachedTo );
		return;
	}

	axpy( &(pActingOn->velocity) , pao.mechanics.attachment_elastic_constant , displacement );

	return;
}
