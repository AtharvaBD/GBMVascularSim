#include "PhysiMeSS_fibre.h"
#include "PhysiMeSS_cell.h"
#include <algorithm>

bool isFibre(PhysiCell::Cell* pCell) 
{
    const auto agentname = std::string(pCell->type_name);
    const auto ecm    = std::string("ecm");
    const auto matrix = std::string("matrix");
    const auto fiber  = std::string("fiber");
    const auto fibre  = std::string("fibre");
    const auto rod    = std::string("rod");

    return (agentname.find(ecm)    != std::string::npos ||
            agentname.find(matrix) != std::string::npos ||
            agentname.find(fiber)  != std::string::npos ||
            agentname.find(fibre)  != std::string::npos ||
            agentname.find(rod)    != std::string::npos);
}

bool isFibre(PhysiCell::Cell_Definition * cellDef)
{
    const auto agentname = std::string(cellDef->name);
    const auto ecm    = std::string("ecm");
    const auto matrix = std::string("matrix");
    const auto fiber  = std::string("fiber");
    const auto fibre  = std::string("fibre");
    const auto rod    = std::string("rod");

    return (agentname.find(ecm)    != std::string::npos ||
            agentname.find(matrix) != std::string::npos ||
            agentname.find(fiber)  != std::string::npos ||
            agentname.find(fibre)  != std::string::npos ||
            agentname.find(rod)    != std::string::npos);
}

std::vector<PhysiCell::Cell_Definition*>* getFibreCellDefinitions() {
    std::vector<PhysiCell::Cell_Definition*>* result = new std::vector<PhysiCell::Cell_Definition*>();
    PhysiCell::Cell_Definition* pCD;

    for (auto& cd_name: PhysiCell::cell_definitions_by_name) {
        if (isFibre(cd_name.second)) {
            result->push_back(cd_name.second);
        }
    }

    return result;
}

PhysiMeSS_Fibre::PhysiMeSS_Fibre() 
{
    fibres_crosslinkers.clear();
    fibres_crosslink_point.clear();

    // length / radius will be set in assign_fibre_orientation()
    X_crosslink_count = 0;
    fail_count = 0;
}

void PhysiMeSS_Fibre::assign_fibre_orientation() 
{ 
    // Enforce total fibre length = 25 µm → half-length = 12.5
    mLength = 12.5;  // microns (half-length)
    mRadius = this->custom_data["fibre_radius"];

    // Initial orientation (may be overridden below)
    this->assign_orientation();

    if (default_microenvironment_options.simulate_2D) {
        // Special case: ECM fibres tangential to circle through their center
        if (this->type_name.find("ecm") != std::string::npos)
        {
            double x = this->position[0];
            double y = this->position[1];

            // tangential direction (-y, x)
            double tx = -y;
            double ty =  x;
            double norm = sqrt(tx*tx + ty*ty);
            if (norm > 0.0)
            {
                tx /= norm;
                ty /= norm;
            }
            this->state.orientation[0] = tx;
            this->state.orientation[1] = ty;
            this->state.orientation[2] = 0.0;
        }
        else if (this->custom_data["anisotropic_fibres"] > 0.5){
            double theta = PhysiCell::NormalRandom(this->custom_data["fibre_angle"],
                                                   this->custom_data["angle_normdist_sd"]);
            this->state.orientation[0] = cos(theta);
            this->state.orientation[1] = sin(theta);
            this->state.orientation[2] = 0.0;
        }
        else{
            this->state.orientation = PhysiCell::UniformOnUnitCircle();
            this->state.orientation[2] = 0.0;
        }
    }
    else {
        this->state.orientation = PhysiCell::UniformOnUnitSphere();
    }

    //###########################################//
    //   this bit a hack for PacMan and maze     //
    //###########################################//
    if (this->type_name == "fibre_vertical") {
        this->state.orientation[0] = 0.0;
        this->state.orientation[1] = 1.0;
        this->state.orientation[2] = 0.0;
    }
    if (this->type_name == "fibre_horizontal") {
        this->state.orientation[0] = 1.0;
        this->state.orientation[1] = 0.0;
        this->state.orientation[2] = 0.0;
    }
    //###########################################// 
}

void PhysiMeSS_Fibre::check_out_of_bounds(std::vector<double>& position)
{
    double Xmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[0]; 
    double Ymin = BioFVM::get_default_microenvironment()->mesh.bounding_box[1]; 
    double Zmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[2]; 

    double Xmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[3]; 
    double Ymax = BioFVM::get_default_microenvironment()->mesh.bounding_box[4]; 
    double Zmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[5]; 
    
    if( default_microenvironment_options.simulate_2D == true )
    {
        Zmin = 0.0; 
        Zmax = 0.0; 
    }
    
    // start and end points of a fibre are calculated from fibre center
    double xs = position[0] - this->mLength * this->state.orientation[0];
    double xe = position[0] + this->mLength * this->state.orientation[0];
    double ys = position[1] - this->mLength * this->state.orientation[1];
    double ye = position[1] + this->mLength * this->state.orientation[1];
    double zs = 0.0;
    double ze = 0.0;
    if (default_microenvironment_options.simulate_2D) {
        // 2D endpoints logging omitted
    }
    else if (!default_microenvironment_options.simulate_2D) {
        zs = position[2] - this->mLength * this->state.orientation[2];
        ze = position[2] + this->mLength * this->state.orientation[2];
    }

    if (this->custom_data["anisotropic_fibres"]) {
        if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
            ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
            fail_count = 10;
        }
    }
    else{
        if (default_microenvironment_options.simulate_2D) {
            while (fail_count < 10) {
                if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                    ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
                    fail_count++;
                    this->state.orientation = PhysiCell::UniformOnUnitCircle();
                    xs = position[0] - mLength * this->state.orientation[0];
                    xe = position[0] + mLength * this->state.orientation[0];
                    ys = position[1] - mLength * this->state.orientation[1];
                    ye = position[1] + mLength * this->state.orientation[1];
                }
                else {
                    break;
                }
            }
        }

        if (!default_microenvironment_options.simulate_2D) {
            while (fail_count < 10) {
                if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                    ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
                    zs < Zmin || ze > Zmax || ze < Xmin || zs > Xmax) {
                    fail_count++;
                    this->state.orientation = PhysiCell::UniformOnUnitSphere();
                    xs = position[0] - mLength * this->state.orientation[0];
                    xe = position[0] + mLength * this->state.orientation[0];
                    ys = position[1] - mLength * this->state.orientation[1];
                    ye = position[1] + mLength * this->state.orientation[1];
                    zs = position[2] - mLength * this->state.orientation[2];
                    ze = position[2] + mLength * this->state.orientation[2];
                }
                else {
                    break;
                }
            }
        }
    }
}

void PhysiMeSS_Fibre::add_potentials_from_cell(PhysiMeSS_Cell* cell) 
{
    // fibres only get pushed or rotated by motile cells
    if (!cell->phenotype.motility.is_motile || X_crosslink_count >= 2) {
        return;
    }

    double distance = 0.0;
    nearest_point_on_fibre(cell->position, displacement);
    for (int index = 0; index < 3; index++) {
        distance += displacement[index] * displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);
    // fibre should only interact with cell if it comes within cell radius plus fibre radius
    double R = phenotype.geometry.radius + mRadius;
    if (distance <= R) {
        std::vector<double> point_of_impact(3, 0.0);
        for (int index = 0; index < 3; index++) {
            point_of_impact[index] = (*cell).position[index] - displacement[index];
        }
        // cell-fibre pushing only if fibre no crosslinks
        if (X_crosslink_count == 0) {
            // fibre pushing turned on
            if (cell->custom_data["fibre_pushing"] > 0.5) {
                static double simple_pressure_scale = 0.027288820670331;
                double temp_r = 0;
                temp_r = -distance;
                temp_r /= R;
                temp_r += 1.0;
                temp_r *= temp_r;
                state.simple_pressure += (temp_r / simple_pressure_scale);

                double effective_repulsion = sqrt(phenotype.mechanics.cell_cell_repulsion_strength *
                                                  (*cell).phenotype.mechanics.cell_cell_repulsion_strength);
                temp_r *= effective_repulsion;

                if (fabs(temp_r) < 1e-16) { return; }
                temp_r /= distance;
                naxpy(&velocity, temp_r, displacement);
            }

            // fibre rotation turned on (2D)
            if (cell->custom_data["fibre_rotation"] > 0.5) {
                std::vector<double> old_orientation(3, 0.0);
                for (int i = 0; i < 2; i++) {
                    old_orientation[i] = state.orientation[i];
                }

                double moment_arm_magnitude = sqrt(
                        point_of_impact[0] * point_of_impact[0] + point_of_impact[1] * point_of_impact[1]);
                double impulse = cell->custom_data["fibre_sticky"]*(*cell).phenotype.motility.migration_speed * moment_arm_magnitude;
                double fibre_length = 2 * mLength;
                double angular_velocity = impulse / (0.5 * fibre_length * fibre_length);
                double angle = angular_velocity;
                state.orientation[0] = old_orientation[0] * cos(angle) - old_orientation[1] * sin(angle);
                state.orientation[1] = old_orientation[0] * sin(angle) + old_orientation[1] * cos(angle);
                normalize(&state.orientation);
            }
        }

        // fibre rotation around other fibre (2D only and fibres intersect at a single point)
        if (cell->custom_data["fibre_rotation"] > 0.5 && X_crosslink_count == 1) {
            double distance_fibre_centre_to_crosslink = 0.0;
            std::vector<double> fibre_centre_to_crosslink(3, 0.0);
            for (int i = 0; i < 2; i++) {
                fibre_centre_to_crosslink[i] = fibres_crosslink_point[i]-position[i];
                distance_fibre_centre_to_crosslink += fibre_centre_to_crosslink[i]*fibre_centre_to_crosslink[i];
            }
            distance_fibre_centre_to_crosslink = sqrt(distance_fibre_centre_to_crosslink);

            std::vector<double> old_orientation(3, 0.0);
            for (int i = 0; i < 2; i++) {
                old_orientation[i] = state.orientation[i];
            }
            double moment_arm_magnitude = sqrt(
                    point_of_impact[0] * point_of_impact[0] + point_of_impact[1] * point_of_impact[1]);
            double impulse = cell->custom_data["fibre_sticky"]*(*cell).phenotype.motility.migration_speed * moment_arm_magnitude;
            double fibre_length = 2 * mLength;
            double angular_velocity = impulse / (0.5 * fibre_length * fibre_length);
            double angle = angular_velocity;
            state.orientation[0] = old_orientation[0] * cos(angle) - old_orientation[1] * sin(angle);
            state.orientation[1] = old_orientation[0] * sin(angle) + old_orientation[1] * cos(angle);
            normalize(&state.orientation);
            position[0] = fibres_crosslink_point[0]-distance_fibre_centre_to_crosslink*state.orientation[0];
            position[1] = fibres_crosslink_point[1]-distance_fibre_centre_to_crosslink*state.orientation[1];
        }
    }

    return;
}

void PhysiMeSS_Fibre::add_potentials_from_fibre(PhysiMeSS_Fibre* other_fibre) 
{
    // placeholder for future fibre–fibre tension modelling
    return;
}

void PhysiMeSS_Fibre::register_fibre_voxels() 
{
    int voxel;
    int voxel_size = this->get_container()->underlying_mesh.dx; // mechanics_voxel_size
    int test = 2.0 * this->mLength / voxel_size; // sample along the fibre

    std::vector<double> fibre_start(3, 0.0);
    std::vector<double> fibre_end(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        fibre_start[i] = this->position[i] - this->mLength * this->state.orientation[i];
        fibre_end[i]   = this->position[i] + this->mLength * this->state.orientation[i];
    }
    // add voxel of fibre end point
    voxel = this->get_container()->underlying_mesh.nearest_voxel_index(fibre_end);
    physimess_voxels.push_back(voxel);
    if (std::find(this->get_container()->agent_grid[voxel].begin(),
                  this->get_container()->agent_grid[voxel].end(),
                  this) == this->get_container()->agent_grid[voxel].end()) {
        this->get_container()->agent_grid[voxel].push_back(this);
    }
    // walk along fibre
    std::vector<double> point_on_fibre(3, 0.0);
    for (unsigned int j = 0; j < test + 1; j++) {
        for (unsigned int i = 0; i < 3; i++) {
            point_on_fibre[i] = fibre_start[i] + j * voxel_size * this->state.orientation[i];
        }
        voxel = this->get_container()->underlying_mesh.nearest_voxel_index(point_on_fibre);
        physimess_voxels.push_back(voxel);
        if (std::find(this->get_container()->agent_grid[voxel].begin(),
                      this->get_container()->agent_grid[voxel].end(),
                      this) == this->get_container()->agent_grid[voxel].end()) {
            this->get_container()->agent_grid[voxel].push_back(this);
        }
    }

    physimess_voxels.sort();
    physimess_voxels.unique();
}

void PhysiMeSS_Fibre::deregister_fibre_voxels() 
{
    int centre_voxel = this->get_container()->underlying_mesh.nearest_voxel_index(this->position);
    for (int voxel: physimess_voxels) {
        if (voxel != centre_voxel) {
            this->get_container()->remove_agent_from_voxel(this, voxel);
        }
    }
}

std::vector<double> PhysiMeSS_Fibre::nearest_point_on_fibre(std::vector<double> point, std::vector<double> &displacement) 
{
    // don't bother if not a fibre
    if (!isFibre(this)) { return displacement; }

    double fibre_length = 2 * this->mLength;
    std::vector<double> fibre_to_agent(3, 0.0);
    double fibre_to_agent_length_squared = 0;
    double fibre_to_agent_dot_fibre_vector = 0;

    for (unsigned int i = 0; i < 3; i++) {
        fibre_to_agent[i] = point[i] - (this->position[i]
                                        - this->mLength * this->state.orientation[i]);
        fibre_to_agent_length_squared += fibre_to_agent[i] * fibre_to_agent[i];
        fibre_to_agent_dot_fibre_vector += fibre_to_agent[i] * fibre_length * this->state.orientation[i];
    }

    if (fibre_to_agent_dot_fibre_vector < 0.) {
        for (int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i];
        }
    }
    else if (fibre_to_agent_dot_fibre_vector > fibre_length * fibre_length) {
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = point[i] - (this->position[i]
                                          + this->mLength * this->state.orientation[i]);
        }
    }
    else {
        double fibre_to_agent_length_cos_alpha_squared =
                fibre_to_agent_dot_fibre_vector * fibre_to_agent_dot_fibre_vector /
                (fibre_length * fibre_length);
        double l = sqrt(fibre_to_agent_length_cos_alpha_squared);
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i] - l * this->state.orientation[i];
        }
    }

    return displacement;
}

void PhysiMeSS_Fibre::check_fibre_crosslinks(PhysiMeSS_Fibre *fibre_neighbor) 
{
    if (this == fibre_neighbor) { return; }

    if (isFibre(this) && isFibre(fibre_neighbor)) {

        std::vector<double> point1(3, 0.0);
        std::vector<double> point2(3, 0.0);
        std::vector<double> point3(3, 0.0);
        std::vector<double> point4(3, 0.0);
        for (int i = 0; i < 3; i++) {
            point1[i] = this->position[i] - mLength * this->state.orientation[i];
            point2[i] = this->position[i] + mLength * this->state.orientation[i];
            point3[i] = fibre_neighbor->position[i] - fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
            point4[i] = fibre_neighbor->position[i] + fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
        }

        std::vector<double> p1_to_p2(3, 0.0);
        std::vector<double> p3_to_p4(3, 0.0);
        std::vector<double> p1_to_p3(3, 0.0);
        std::vector<double> centre_to_centre(3, 0.0);
        for (int i = 0; i < 3; i++) {
            p1_to_p2[i] = point2[i] - point1[i];
            p3_to_p4[i] = point4[i] - point3[i];
            p1_to_p3[i] = point3[i] - point1[i];
            centre_to_centre[i] = fibre_neighbor->position[i] - this->position[i];
        }

        double co_radius = this->mRadius + fibre_neighbor->mRadius;
        double co_length = this->mLength + fibre_neighbor->mLength;
        std::vector<double> zero(3, 0.0);
        double distance = PhysiCell::dist(zero, centre_to_centre);
        normalize(&centre_to_centre);

        std::vector<double> FCP = cross_product(p1_to_p2, p3_to_p4);
        if (dot_product(FCP,FCP) == 0 &&
            (centre_to_centre == this->state.orientation ||
             centre_to_centre == -1.0 * this->state.orientation) &&
            distance <= co_length) {
            if (std::find(this->fibres_crosslinkers.begin(), this->fibres_crosslinkers.end(), (fibre_neighbor)) ==
                this->fibres_crosslinkers.end()) {
                this->fibres_crosslinkers.push_back(fibre_neighbor);
            }
            this->fibres_crosslink_point = this->position + this->mLength * centre_to_centre;
        }

        std::vector<double> displacement(3, 0.0);
        fibre_neighbor->nearest_point_on_fibre(point1, displacement);
        double test_point1 = PhysiCell::dist(zero, displacement);
        fibre_neighbor->nearest_point_on_fibre(point2, displacement);
        double test_point2 = PhysiCell::dist(zero, displacement);
        this->nearest_point_on_fibre(point3, displacement);
        double test_point3 = PhysiCell::dist(zero, displacement);
        this->nearest_point_on_fibre(point4, displacement);
        double test_point4 = PhysiCell::dist(zero, displacement);
        if (std::abs(test_point1) <= co_radius ||
            std::abs(test_point2) <= co_radius ||
            std::abs(test_point3) <= co_radius ||
            std::abs(test_point4) <= co_radius &&
            centre_to_centre != this->state.orientation &&
            centre_to_centre != -1.0 * this->state.orientation) {
            if (std::find(this->fibres_crosslinkers.begin(), this->fibres_crosslinkers.end(), (fibre_neighbor)) ==
                this->fibres_crosslinkers.end()) {
                this->fibres_crosslinkers.push_back(fibre_neighbor);
            }
            this->fibres_crosslink_point = point1;
        }

        double test2_tolerance = co_radius;
        double test2 = dot_product(p1_to_p3, FCP);
        if (std::abs(test2) < test2_tolerance) {
            double a = dot_product(p1_to_p2, p1_to_p3) / dot_product(p1_to_p2, p1_to_p2);
            double b = dot_product(p1_to_p2, p3_to_p4) / dot_product(p1_to_p2, p1_to_p2);
            std::vector<double> c(3, 0.0);
            std::vector<double> n(3, 0.0);
            for (int i = 0; i < 3; i++) {
                c[i] = b * p1_to_p2[i] - p3_to_p4[i];
                n[i] = p1_to_p3[i] - a * p1_to_p2[i];
            }
            double t_neighbor = dot_product(c, n) / dot_product(c, c);
            double t_this = a + b * t_neighbor;
            std::vector<double> crosslink_point(3, 0.0);
            for (int i = 0; i < 2; i++) {
                crosslink_point[i] = point1[i] + t_this * p1_to_p2[i];
            }
            double tolerance = co_radius/co_length;
            double lower_bound = 0.0 - tolerance;
            double upper_bound = 1.0 + tolerance;
            if (lower_bound <= t_neighbor && t_neighbor <= upper_bound &&
                lower_bound <= t_this && t_this <= upper_bound) {
                if (std::find(this->fibres_crosslinkers.begin(), this->fibres_crosslinkers.end(), (fibre_neighbor)) ==
                    this->fibres_crosslinkers.end()) {
                    this->fibres_crosslinkers.push_back(fibre_neighbor);
                }
                this->fibres_crosslink_point = crosslink_point;
            }
        }
    } else { return; }
}

void PhysiMeSS_Fibre::add_crosslinks() 
{
    for (auto* neighbor : physimess_neighbors)
    {
        if (isFibre(neighbor)) {
            this->check_fibre_crosslinks(static_cast<PhysiMeSS_Fibre*>(neighbor));
        }
    }
}
