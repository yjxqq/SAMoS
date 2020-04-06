/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file pair_vertex_particle_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Implementation of PairVertexParticlePotential class
 */

#include "my_pair_vertex_particle_potential.hpp"
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_expression.hpp>
//#include <boost/numeric/ublas/traits/c_array.hpp>
//#include <boost/numeric/ublas/assignment.hpp>

using namespace boost::numeric::ublas;

static int prev(int i, int N)
{
    if (i == 0) return N - 1;
    return i - 1;
}

static int next(int i, int N)
{
    if (i == N - 1) return 0;
    return i + 1;
}

//! \param dt time step sent by the integrator 
void MyPairVertexParticlePotential::compute(double dt)
{
    int N = m_system->size();
    double K = m_K;
    double gamma = m_gamma;
    double lambda = m_lambda;
    double alpha = 1.0;  // phase in factor
    double pot_eng = 0.0;

    // My changes
    bool consider_edge_myosin_activity = false;
    bool consider_perim_myosin_activity = false;
    bool consider_cell_substrate_adhesion = true;

    double H_alpha = 2.0;
    double H_base = 0.0;
    double H_n = 8.0;
    double H_Ps = 6.73;
    double H_KL = 1.0;
    double H_Ls = 6.73 / 6.0;

    if (m_mesh_update_steps > 0) if (m_system->get_step() % m_mesh_update_steps == 0) m_nlist->build_mesh();

    Mesh& mesh = m_system->get_mesh();

    if (m_system->compute_per_particle_energy())
    {
        for (int i = 0; i < N; i++)
        {
            Particle& p = m_system->get_particle(i);
            p.set_pot_energy("my_vp", 0.0);
        }
    }

    m_potential_energy = 0.0;
    for (int i = 0; i < N; i++)
    {
        Particle& pi = m_system->get_particle(i);
        Vertex& vi = mesh.get_vertices()[i];
        Vector3d Nvec = Vector3d(pi.Nx, pi.Ny, pi.Nz);
        if (m_phase_in) alpha = m_val->get_val(static_cast<int>(pi.age / dt));
        // ##############First handle the cell itself ######################
        if (pi.in_tissue && (m_include_boundary || !vi.boundary))
        {
            if (m_has_part_params)
            {
                K = m_particle_params[vi.type - 1].K;
                gamma = m_particle_params[vi.type - 1].gamma;
                lambda = m_particle_params[vi.type - 1].lambda;
            }
            double dA = 0.0;
            dA = (vi.area - pi.A0);
            double area_term = 0.5 * K * dA;
            double perim_term = gamma * vi.perim;
            double perim_term_edge_myosin = 0.0;
            double perim_term_perim_myosin = gamma * vi.perim
                    * (H_base + H_alpha * (pow(vi.perim / H_Ps, H_n)) / (H_KL + pow(vi.perim / H_Ps, H_n)));
            pot_eng = 0.5 * (K * dA * dA + gamma * vi.perim * vi.perim);
            Vector3d area_vec(0.0, 0.0, 0.0);
            Vector3d perim_vec(0.0, 0.0, 0.0);
            Vector3d con_vec(0.0, 0.0, 0.0);
            Vector3d sub_adh_vec(0.0, 0.0, 0.0);
            for (int f = 0; f < vi.n_faces; f++)
            {
                int fid = vi.dual[f];
                int fid_m = vi.dual[prev(f, vi.n_faces)];    // why using vi.dual? what is dual?
                int fid_p = vi.dual[next(f, vi.n_faces)];    // dual should be the dual vertices of a cell.

                Face& f_nu_m = mesh.get_faces()[fid_m];
                Face& f_nu = mesh.get_faces()[fid];
                Face& f_nu_p = mesh.get_faces()[fid_p];

                Vector3d r_nu_m = f_nu_m.rc;
                Vector3d r_nu = f_nu.rc;
                Vector3d r_nu_p = f_nu_p.rc;

                Vector3d cross_prod_1(0.0, 0.0, 0.0);
                if (!(f_nu.is_hole || f_nu_p.is_hole))
                { /*std::cout << " is_hole ";*/
                    cross_prod_1 = (cross(r_nu_p, Nvec)) * f_nu.get_jacobian(i);
                }
                if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_1 = cross_prod_1
                        - (cross(r_nu_m, Nvec)) * f_nu.get_jacobian(i);
                area_vec = area_vec + cross_prod_1;

                // My changes
                Vector3d cross_prod_2(0.0, 0.0, 0.0);
                double mij_p = H_base
                        + H_alpha * pow((r_nu_p - r_nu).len() / H_Ls, H_n)
                                / (H_KL + pow((r_nu_p - r_nu).len() / H_Ls, H_n));

                if (!(f_nu_m.is_hole || f_nu.is_hole))
                {
                    if (consider_edge_myosin_activity)
                    {
                        double mij_m = H_base
                                + H_alpha * pow((r_nu - r_nu_m).len() / H_Ls, H_n)
                                        / (H_KL + pow((r_nu - r_nu_m).len() / H_Ls, H_n));
                        cross_prod_2 = (sqrt(mij_m) * (r_nu - r_nu_m).unit()) * f_nu.get_jacobian(i);
                    }
                    else cross_prod_2 = ((r_nu - r_nu_m).unit()) * f_nu.get_jacobian(i);
                }
                if (!(f_nu_p.is_hole || f_nu.is_hole))
                {
                    if (consider_edge_myosin_activity)
                    {
                        perim_term_edge_myosin += gamma * sqrt(mij_p) * (r_nu_p - r_nu).len();
                        double mij_p = H_base
                                + H_alpha * pow((r_nu_p - r_nu).len() / H_Ls, H_n)
                                        / (H_KL + pow((r_nu_p - r_nu).len() / H_Ls, H_n));
                        cross_prod_2 -= (sqrt(mij_p) * (r_nu_p - r_nu).unit()) * f_nu.get_jacobian(i);
                    }
                    else cross_prod_2 -= ((r_nu_p - r_nu).unit()) * f_nu.get_jacobian(i);
                }
                perim_vec = perim_vec + cross_prod_2;

                // Adding force contributions from the dual edges associated with face f and face prev(f,vi.n_faces)
                if (m_has_pair_params)
                {
                    Vertex& vn = mesh.get_vertices()[vi.dual_neighbour_map[f]];
                    lambda = m_pair_params[vi.type - 1][vn.type - 1].lambda;
                }
                Vector3d cross_prod_3(0.0, 0.0, 0.0);
                if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_3 = lambda
                        * (((r_nu - r_nu_m).unit()) * f_nu.get_jacobian(i));
                con_vec = con_vec + cross_prod_3;

                // Adding force contributions from the dual edges associated with face f and face next(f,vi.n_faces)
                if (m_has_pair_params)
                {
                    Vertex& vn = mesh.get_vertices()[vi.dual_neighbour_map[next(f, vi.n_faces)]];
                    lambda = m_pair_params[vi.type - 1][vn.type - 1].lambda;
                }
                Vector3d cross_prod_4(0.0, 0.0, 0.0);
                if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_4 = lambda
                        * (((r_nu_p - r_nu).unit()) * f_nu.get_jacobian(i));
                con_vec = con_vec - cross_prod_4;

                pot_eng += lambda * (r_nu - r_nu_m).len();

                /*--------------------------------------Cell Area adhesion contribution STRAT---------------------------------------------*/
                /*---------------------------------------------------148----678------------------------------------------------------------*/
                /*----------------------------------------SETTING PARAMETERS FIRST--------------------------------------------------*/

                double stripe_width = 5;
                double stripe_distance = 30.0;

                //double initial_xcoord = -45.0;
                double initial_xcoord = 24.5;
                double symmetry_line_xcoord = -0.0;

                //bool IsSymmetric = true;
                bool IsSymmetric = false;

                double adhesive_area_parameter = -1.0;
                //double adhesive_area_parameter = 0.0;

                //double small_change = 0.02;
                //double preferred_sample_dis = 0.004;

                //double small_change = 0.01;
                //double preferred_sample_dis = 0.002;

                double small_change = 0.1;
                double preferred_sample_dis = 0.02;

                //double small_change = 0.2;
                //double preferred_sample_dis = 0.04;

                bool adhesive_area_parameter_can_change = false;
                //bool adhesive_area_parameter_can_change = true;

                /*--------------------------   ----------- PARAMETERS SETTING DONE--------------------------------------------------*/
                if (consider_cell_substrate_adhesion)
                {
                    bool is_symm_counterpart = false;
                    c_vector<double, 2> area_adhesion_contribution = zero_vector<double>(2);
                    c_vector<double, 2> adhesive_area_gradient = zero_vector<double>(2);
                    c_vector<double, 2> weighted_adhesive_area_gradient = zero_vector<double>(2);

                    int index = -1;
                    for (int k = 0; k < 3; k++)
                    {
                        if (f_nu.vertices[k] == i)
                        {
                            index = k;
                            break;
                        }
                    }
                    if (index == -1) std::cerr << "err00" << '\n';
                    Vertex& vm = mesh.get_vertices()[f_nu.vertices[(index + 2) % 3]];
                    Vertex& vn = mesh.get_vertices()[f_nu.vertices[(index + 1) % 3]];

                    int n_real_cell_attached = 0;
                    if ((!vm.boundary) && (!vn.boundary)) n_real_cell_attached = 3;
                    else if (vm.boundary && vn.boundary) n_real_cell_attached = 1;
                    else
                    {
                        n_real_cell_attached = 2;
                    }
                    if (n_real_cell_attached == 0) std::cerr
                            << "Error0: face center does not attached to any real cell";

                    if (n_real_cell_attached != 3)
                    {
                        if (n_real_cell_attached == 2)
                        {
                            /*                if (vm.boundary)
                             {
                             r_nu_p = r_nu_m;
                             int g = -1;
                             for (int fn = 0; fn < vn.n_faces; fn++)
                             {
                             if (vn.dual[fn]==fid)
                             {
                             g= fn;
                             break;
                             }
                             }
                             r_nu_m = (mesh.get_faces()[vn.dual[next(g,vn.n_faces)]]).rc;

                             if (g == -1)
                             {
                             std::cerr << "Err2" << r_nu_m.x <<'\t' << r_nu_m.y << '\t' << r_nu.x <<'\t' << r_nu.y << '\t' << r_nu_p.x <<'\t' << r_nu_p.y << '\t' << '\n';
                             std::cout << "Situation_second " << '\n';
                             std::cout << "f_nu.vertices: " << fid << ' ' << mesh.get_vertices()[f_nu.vertices[0]] << ' ' << mesh.get_vertices()[f_nu.vertices[1]] << ' ' << mesh.get_vertices()[f_nu.vertices[2]] << '\n';
                             std::cout << "vertices now using: " << vi.dual_neighbour_map[f] << ' ' << i <<' ' << vi.dual_neighbour_map[next(f,vi.n_faces)] << '\n';
                             for (int f = 0; f < vi.n_faces; f++)
                             {
                             std::cout << vi.dual[f] << '\t';
                             }
                             std::cout << '\n';

                             for  (int fn = 0; fn < vn.n_faces; fn++)
                             {
                             std::cout << vn.dual[fn] << '\t';
                             }
                             std::cout << '\n';
                             for (int fm = 0; fm < vm.n_faces; fm++)
                             {
                             std::cout << vm.dual[fm] << '\t';
                             }
                             std::cout << '\n';
                             }
                             }
                             else
                             {
                             r_nu_m = r_nu_p;
                             int g = -1;
                             for (int fm = 0; fm < vm.n_faces; fm++)
                             {
                             if (vm.dual[fm]==fid)
                             {
                             g= fm;
                             break;
                             }
                             }
                             r_nu_p = (mesh.get_faces()[vm.dual[prev(g,vm.n_faces)]]).rc;

                             if (g == -1)
                             std::cerr << "Err1";
                             }*/
                            if (vm.boundary)
                            {
                                int g = -1;
                                for (int fm = 0; fm < vm.n_faces; fm++)
                                {
                                    if (vm.faces[fm] == fid)
                                    {
                                        g = fm;
                                        break;
                                    }
                                }
                                r_nu_m = mesh.get_faces()[vm.faces[prev(g, vm.n_faces)]].rc;
                                r_nu_p = mesh.get_faces()[vm.faces[next(g, vm.n_faces)]].rc;
                                if (g == -1) std::cerr << "err1" << '\n';
                            }
                            else
                            {
                                int g = -1;
                                for (int fn = 0; fn < vn.n_faces; fn++)
                                {
                                    if (vn.faces[fn] == fid)
                                    {
                                        g = fn;
                                        break;
                                    }
                                }
                                r_nu_m = mesh.get_faces()[vn.faces[prev(g, vn.n_faces)]].rc;
                                r_nu_p = mesh.get_faces()[vn.faces[next(g, vn.n_faces)]].rc;
                                if (g == -1) std::cerr << "err2" << '\n';
                            }
                        } // End of altering the typical 3 points

                        double triangle_box_bottom = 0.0;
                        double triangle_box_top = 0.0;
                        double triangle_box_left = 0.0;
                        double triangle_box_right = 0.0;
                        double sample_area_bottom = 0.0;
                        double sample_area_top = 0.0;
                        double sample_area_left = 0.0;
                        double sample_area_right = 0.0;
                        double sample_area = 0.0;
                        double sample_area_bottom_x = 0.0;
                        double sample_area_top_x = 0.0;
                        double sample_area_left_x = 0.0;
                        double sample_area_right_x = 0.0;
                        double sample_area_x = 0.0;
                        double sample_area_bottom_y = 0.0;
                        double sample_area_top_y = 0.0;
                        double sample_area_left_y = 0.0;
                        double sample_area_right_y = 0.0;
                        double sample_area_y = 0.0;
                        unsigned sample_num = 0;
                        unsigned sample_num_x = 0;
                        unsigned sample_num_y = 0;
                        unsigned num_across = 0;
                        unsigned num_up = 0;
                        unsigned num_across_x = 0;
                        unsigned num_up_x = 0;
                        unsigned num_across_y = 0;
                        unsigned num_up_y = 0;

                        c_vector<double, 2> centroid = zero_vector<double>(2);
                        centroid[0] = (r_nu_m.x + r_nu.x + r_nu_p.x) / 3;
                        centroid[1] = (r_nu_m.y + r_nu.y + r_nu_p.y) / 3;
                        int floornum = floor(centroid[1] / stripe_distance);
                        //unsigned stripe_num = centroid[1]<(floornum+0.5)*stripe_distance ? floornum:floornum+1;
                        int stripe_num = centroid[1] < (floornum + 0.5) * stripe_distance ? floornum : floornum + 1;
                        double stripe_location = stripe_distance * stripe_num;
                        double stripe_bottom = stripe_location - stripe_width / 2;
                        double stripe_top = stripe_location + stripe_width / 2;

                        // Preparation for the original triangle and its triangle box and adhesive_area box
                        bool has_adhesive_area = true;
                        triangle_box_bottom = std::min(std::min(r_nu_m.y, r_nu.y), r_nu_p.y);
                        triangle_box_top = std::max(std::max(r_nu_m.y, r_nu.y), r_nu_p.y);
                        triangle_box_left = std::min(std::min(r_nu_m.x, r_nu.x), r_nu_p.x);
                        triangle_box_right = std::max(std::max(r_nu_m.x, r_nu.x), r_nu_p.x);
                        // Changes for symmetry:
                        if (IsSymmetric == true)
                        {
                            if ((triangle_box_left + triangle_box_right) / 2.0 < symmetry_line_xcoord)
                            {
                                is_symm_counterpart = true;
                                double new_triangle_box_left = symmetry_line_xcoord * 2.0 - triangle_box_right;
                                double new_triangle_box_right = symmetry_line_xcoord * 2.0 - triangle_box_left;
                                triangle_box_left = new_triangle_box_left;
                                triangle_box_right = new_triangle_box_right;
                            }
                        }

                        if (triangle_box_right < initial_xcoord) has_adhesive_area = false;
                        else
                        {
                            sample_area_right = triangle_box_right;
                            sample_area_left = triangle_box_left > initial_xcoord ? triangle_box_left : initial_xcoord;
                        }
                        // Changes for symmetry:
                        if (IsSymmetric == true)
                        {
                            if (is_symm_counterpart == true)
                            {
                                if (has_adhesive_area == true)
                                {
                                    double new_sample_area_left = symmetry_line_xcoord * 2.0 - sample_area_right;
                                    double new_sample_area_right = symmetry_line_xcoord * 2.0 - sample_area_left;
                                    sample_area_left = new_sample_area_left;
                                    sample_area_right = new_sample_area_right;
                                }
                            }
                        }

                        if (triangle_box_top < stripe_bottom) //previous error: if (triangle_box_top<triangle_box_bottom)
                        has_adhesive_area = false;
                        else if (triangle_box_bottom < stripe_bottom)
                        {
                            sample_area_bottom = stripe_bottom;
                            if (triangle_box_top < stripe_top) sample_area_top = triangle_box_top;
                            else sample_area_top = stripe_top;
                        }
                        else if (triangle_box_top < stripe_top)
                        {
                            sample_area_bottom = triangle_box_bottom;
                            sample_area_top = triangle_box_top;
                        }
                        else if (triangle_box_bottom < stripe_top)
                        {
                            sample_area_bottom = triangle_box_bottom;
                            sample_area_top = stripe_top;
                        }
                        else has_adhesive_area = false;
                        if (has_adhesive_area == true)
                        {
                            num_across = round((sample_area_right - sample_area_left) / preferred_sample_dis);
                            num_up = round((sample_area_top - sample_area_bottom) / preferred_sample_dis);
                            sample_num = num_across * num_up;
                            sample_area = (sample_area_top - sample_area_bottom)
                                    * (sample_area_right - sample_area_left);
                        }

                        // Preparation for the triangle with small change in x and its triangle box and adhesive_area box
                        bool has_adhesive_area_after_small_change_x = true;
                        r_nu.x += small_change;
                        triangle_box_bottom = std::min(std::min(r_nu_m.y, r_nu.y), r_nu_p.y);
                        triangle_box_top = std::max(std::max(r_nu_m.y, r_nu.y), r_nu_p.y);
                        triangle_box_left = std::min(std::min(r_nu_m.x, r_nu.x), r_nu_p.x);
                        triangle_box_right = std::max(std::max(r_nu_m.x, r_nu.x), r_nu_p.x);
                        r_nu.x -= small_change;
                        // Changes for symmetry:
                        if (IsSymmetric == true)
                        {
                            if (is_symm_counterpart == true)
                            {
                                double new_triangle_box_left = symmetry_line_xcoord * 2.0 - triangle_box_right;
                                double new_triangle_box_right = symmetry_line_xcoord * 2.0 - triangle_box_left;
                                triangle_box_left = new_triangle_box_left;
                                triangle_box_right = new_triangle_box_right;
                            }
                        }

                        if (triangle_box_right < initial_xcoord) has_adhesive_area_after_small_change_x = false;
                        else
                        {
                            sample_area_right_x = triangle_box_right;
                            sample_area_left_x =
                                    triangle_box_left > initial_xcoord ? triangle_box_left : initial_xcoord;
                        }
                        // Changes for symmetry:
                        if (IsSymmetric == true)
                        {
                            if (is_symm_counterpart == true)
                            {
                                if (has_adhesive_area_after_small_change_x == true)
                                {
                                    double new_sample_area_left_x = symmetry_line_xcoord * 2.0 - sample_area_right_x;
                                    double new_sample_area_right_x = symmetry_line_xcoord * 2.0 - sample_area_left_x;
                                    sample_area_left_x = new_sample_area_left_x;
                                    sample_area_right_x = new_sample_area_right_x;
                                }
                            }
                        }

                        if (triangle_box_top < stripe_bottom) has_adhesive_area_after_small_change_x = false;
                        else if (triangle_box_bottom < stripe_bottom)
                        {
                            sample_area_bottom_x = stripe_bottom;
                            if (triangle_box_top < stripe_top) sample_area_top_x = triangle_box_top;
                            else sample_area_top_x = stripe_top;
                        }
                        else if (triangle_box_top < stripe_top)
                        {
                            sample_area_bottom_x = triangle_box_bottom;
                            sample_area_top_x = triangle_box_top;
                        }
                        else if (triangle_box_bottom < stripe_top)
                        {
                            sample_area_bottom_x = triangle_box_bottom;
                            sample_area_top_x = stripe_top;
                        }
                        else has_adhesive_area_after_small_change_x = false;
                        if (has_adhesive_area_after_small_change_x)
                        {
                            num_across_x = round((sample_area_right_x - sample_area_left_x) / preferred_sample_dis);
                            num_up_x = round((sample_area_top_x - sample_area_bottom_x) / preferred_sample_dis);
                            sample_num_x = num_across_x * num_up_x;
                            sample_area_x = (sample_area_top_x - sample_area_bottom_x)
                                    * (sample_area_right_x - sample_area_left_x);
                        }

                        // Preparation for the  triangle with small change in y and its triangle box and adhesive_area box
                        bool has_adhesive_area_after_small_change_y = true;
                        r_nu.y += small_change;
                        triangle_box_bottom = std::min(std::min(r_nu_m.y, r_nu.y), r_nu_p.y);
                        triangle_box_top = std::max(std::max(r_nu_m.y, r_nu.y), r_nu_p.y);
                        triangle_box_left = std::min(std::min(r_nu_m.x, r_nu.x), r_nu_p.x);
                        triangle_box_right = std::max(std::max(r_nu_m.x, r_nu.x), r_nu_p.x);
                        r_nu.y -= small_change;
                        // Changes for symmetry:
                        if (IsSymmetric == true)
                        {
                            if (is_symm_counterpart == true)
                            {
                                double new_triangle_box_left = symmetry_line_xcoord * 2.0 - triangle_box_right;
                                double new_triangle_box_right = symmetry_line_xcoord * 2.0 - triangle_box_left;
                                triangle_box_left = new_triangle_box_left;
                                triangle_box_right = new_triangle_box_right;
                            }
                        }

                        if (triangle_box_right < initial_xcoord) has_adhesive_area_after_small_change_y = false;
                        else
                        {
                            sample_area_right_y = triangle_box_right;
                            sample_area_left_y =
                                    triangle_box_left > initial_xcoord ? triangle_box_left : initial_xcoord;
                        }
                        // Changes for symmetry:
                        if (IsSymmetric == true)
                        {
                            if (is_symm_counterpart == true)
                            {
                                if (has_adhesive_area_after_small_change_y == true)
                                {
                                    double new_sample_area_left_y = symmetry_line_xcoord * 2.0 - sample_area_right_y;
                                    double new_sample_area_right_y = symmetry_line_xcoord * 2.0 - sample_area_left_y;
                                    sample_area_left_y = new_sample_area_left_y;
                                    sample_area_right_y = new_sample_area_right_y;
                                }
                            }
                        }

                        if (triangle_box_top < stripe_bottom) has_adhesive_area_after_small_change_y = false;
                        else if (triangle_box_bottom < stripe_bottom)
                        {
                            sample_area_bottom_y = stripe_bottom;
                            if (triangle_box_top < stripe_top) sample_area_top_y = triangle_box_top;
                            else sample_area_top_y = stripe_top;
                        }
                        else if (triangle_box_top < stripe_top)
                        {
                            sample_area_bottom_y = triangle_box_bottom;
                            sample_area_top_y = triangle_box_top;
                        }
                        else if (triangle_box_bottom < stripe_top)
                        {
                            sample_area_bottom_y = triangle_box_bottom;
                            sample_area_top_y = stripe_top;
                        }
                        else has_adhesive_area_after_small_change_y = false;
                        if (has_adhesive_area_after_small_change_y)
                        {
                            num_across_y = round((sample_area_right_y - sample_area_left_y) / preferred_sample_dis);
                            num_up_y = round((sample_area_top_y - sample_area_bottom_y) / preferred_sample_dis);
                            sample_num_y = num_across_y * num_up_y;
                            sample_area_y = (sample_area_top_y - sample_area_bottom_y)
                                    * (sample_area_right_y - sample_area_left_y);
                        }

                        // IF HAS ADHESIVE AREA THEN CALCULATE THE GRADIENT;
                        if (has_adhesive_area || has_adhesive_area_after_small_change_x
                                || has_adhesive_area_after_small_change_y)
                        {
                            c_vector<c_vector<double, 2>, 3> points;
                            points[0][0] = r_nu_m.x;
                            points[0][1] = r_nu_m.y;
                            points[1][0] = r_nu.x;
                            points[1][1] = r_nu.y;
                            points[2][0] = r_nu_p.x;
                            points[2][1] = r_nu_p.y;
                            // Calculate initial adhesive area
                            double weighted_adhesive_sample_num = 0.0;
                            double adhesive_area = 0.0;
                            double weighted_adhesive_area;
                            if (has_adhesive_area && (sample_num > 0))
                            {
                                unsigned adhesive_sample_num = 0;
                                for (unsigned i = 0; i < sample_num; i++)
                                {
                                    bool point_is_outside = false;
                                    c_vector<double, 2> point; // sample point
                                    double x_coord = sample_area_left
                                            + (sample_area_right - sample_area_left) / num_across
                                                    * (i % num_across + 0.5);
                                    double y_coord = sample_area_bottom
                                            + (sample_area_top - sample_area_bottom) / num_up * (i / num_across + 0.5);
                                    point[0] = x_coord;
                                    point[1] = y_coord;

                                    for (unsigned j = 0; j < 3; j++)
                                    {
                                        c_vector<double, 2> vec1 = point - points[j]; //from sample point to point nu
                                        c_vector<double, 2> vec2 = points[(j + 1) % 3] - points[j]; // from point nu to point nu+1

                                        if ((vec1[0] * vec2[1] - vec1[1] * vec2[0]) > 0.0)
                                        {
                                            point_is_outside = true;
                                            break;
                                        }
                                    }
                                    if (point_is_outside == false)
                                    {
                                        adhesive_sample_num += 1;
                                        weighted_adhesive_sample_num += 1.0 + 0.1 * (x_coord - initial_xcoord);
                                    }
                                }
                                adhesive_area = adhesive_sample_num / double(sample_num)
                                        * (sample_area_top - sample_area_bottom)
                                        * (sample_area_right - sample_area_left);
                                weighted_adhesive_area = weighted_adhesive_sample_num / double(sample_num)
                                        * (sample_area_top - sample_area_bottom)
                                        * (sample_area_right - sample_area_left);
                            }

                            // Calculate adhesive area gradient with face center nu;
                            for (unsigned j = 0; j < 2; j++)
                            {
                                double new_adhesive_area = 0.0;
                                double new_weighted_adhesive_area = 0.0;
                                c_vector<c_vector<double, 2>, 3> new_points;
                                new_points = points;
                                c_vector<double, 2> small_vector;
                                small_vector[0] = small_change * (j == 0 ? 1.0 : 0.0);
                                small_vector[1] = small_change * (j == 0 ? 0.0 : 1.0);
                                new_points[1] += small_vector;

                                if ((j == 0 && has_adhesive_area_after_small_change_x && (sample_num_x > 0))
                                        || (j == 1 && has_adhesive_area_after_small_change_y) && (sample_num_y > 0))
                                {
                                    unsigned adhesive_sample_num = 0;
                                    double weighted_adhesive_sample_num = 0.0;
                                    sample_area_left = (j == 0 ? sample_area_left_x : sample_area_left_y);
                                    sample_area_right = (j == 0 ? sample_area_right_x : sample_area_right_y);
                                    sample_area_bottom = (j == 0 ? sample_area_bottom_x : sample_area_bottom_y);
                                    sample_area_top = (j == 0 ? sample_area_top_x : sample_area_top_y);
                                    sample_num = (j == 0 ? sample_num_x : sample_num_y);
                                    num_across = (j == 0 ? num_across_x : num_across_y);
                                    num_up = (j == 0 ? num_up_x : num_up_y);

                                    for (unsigned i = 0; i < sample_num; i++)
                                    {
                                        bool point_is_outside = false;
                                        double x_coord = sample_area_left
                                                + (sample_area_right - sample_area_left) / num_across
                                                        * (i % num_across + 0.5);
                                        double y_coord = sample_area_bottom
                                                + (sample_area_top - sample_area_bottom) / num_up
                                                        * (i / num_across + 0.5);
                                        c_vector<double, 2> point;
                                        point[0] = x_coord;
                                        point[1] = y_coord;
                                        for (unsigned j = 0; j < 3; j++)
                                        {
                                            c_vector<double, 2> vec1 = point - new_points[j];
                                            c_vector<double, 2> vec2 = new_points[(j + 1) % 3] - new_points[j];

                                            if ((vec1[0] * vec2[1] - vec1[1] * vec2[0]) > 0.0)
                                            {
                                                point_is_outside = true;
                                                break;
                                            }
                                        }
                                        if (point_is_outside == false)
                                        {
                                            adhesive_sample_num += 1;
                                            weighted_adhesive_sample_num += 1.0 + 0.1 * (x_coord - initial_xcoord);
                                        }
                                    } // Get adhesive_sample_num

                                    /*
                                     if (adhesive_sample_num!= 0)
                                     std::cout << "good news! " << "adhesive_sample_num!=0 " << adhesive_sample_num <<'\n';
                                     */

                                    new_adhesive_area = adhesive_sample_num / double(sample_num)
                                            * (sample_area_top - sample_area_bottom)
                                            * (sample_area_right - sample_area_left);
                                    new_weighted_adhesive_area = weighted_adhesive_sample_num / double(sample_num)
                                            * (sample_area_top - sample_area_bottom)
                                            * (sample_area_right - sample_area_left);
                                } // end of "if ((j==0 && has_adhesive_a..."; end of claculation of new adhesive area
                                adhesive_area_gradient += (new_adhesive_area - adhesive_area) / small_change
                                        * (small_vector / small_change);
                                weighted_adhesive_area_gradient += (new_weighted_adhesive_area - weighted_adhesive_area)
                                        / small_change * (small_vector / small_change);
                                //std::cout << "new_adhesive_area " << new_adhesive_area << " adhesive_area " << adhesive_area << " adhesive_area_gradient " << adhesive_area_gradient[0] << ' ' << adhesive_area_gradient[1] << '\n';
                            } // end of "for (unsigned j = 0; j<2; j++)"

                        } // end of statement 'if (has_adhesive_area || has_...)' HERE WE GET adhesive_area_gradient

                        if (adhesive_area_parameter_can_change == false) area_adhesion_contribution =
                                -adhesive_area_parameter * adhesive_area_gradient;
                        else area_adhesion_contribution = -adhesive_area_parameter * weighted_adhesive_area_gradient;
                        if (n_real_cell_attached == 2) area_adhesion_contribution = -area_adhesion_contribution;
                    } // END OF "if (n_real_cell_attached!=3)"; END OF Calculation of area_adhesion_contribution from face centre nu done!
                    Vector3d cross_pro_s(0.0, 0.0, 0.0);
                    cross_pro_s.x = area_adhesion_contribution[0];
                    cross_pro_s.y = area_adhesion_contribution[1];
                    cross_pro_s = cross_pro_s * f_nu.get_jacobian(i);
                    sub_adh_vec = sub_adh_vec + cross_pro_s; //pay attention to the +/- sign
                    //END OF MY CALCULATIONS: substrate adhesion force rising from face f
                }
                /*-----------------------------------------------------------------------------------------------------------------------*/
                /*--------------------------------------Cell Area adhesion contribution END---------------------------------------------*/

            } // end of "for (int f = 0; f < vi.n_faces; f++)"

            // area term
            double area_fact = -alpha * area_term;
            pi.fx += area_fact * area_vec.x;
            pi.fy += area_fact * area_vec.y;
            pi.fz += area_fact * area_vec.z;
            // perimeter term
            // My changes
            double perim_fact = -alpha
                    * (consider_edge_myosin_activity ? perim_term_edge_myosin :
                            (consider_perim_myosin_activity ? perim_term_perim_myosin : perim_term));
            pi.fx += perim_fact * perim_vec.x;
            pi.fy += perim_fact * perim_vec.y;
            pi.fz += perim_fact * perim_vec.z;
            // Add contractile term
            pi.fx -= alpha * con_vec.x;
            pi.fy -= alpha * con_vec.y;
            pi.fz -= alpha * con_vec.z;
            // Add substrate adhesion term
            pi.fx += sub_adh_vec.x;
            pi.fy += sub_adh_vec.y;
            pi.fz += sub_adh_vec.z;
            /*
             std::cout << "other forces in x " << area_fact*area_vec.x << ' ' << perim_fact*perim_vec.x << ' ' << '\n';
             if (sub_adh_vec.x != 0.0)
             std::cout << "substrate adhesion force!=0 " << sub_adh_vec.x <<' ' << sub_adh_vec.y << '\n' << std::endl;
             else
             std::cout << "substrate adhesion force=0 " << '\n' << std::endl;
             */

            // add potential energy
            m_potential_energy += pot_eng;
        }
        // ##############Now check neighbours#############
        for (int j = 0; j < vi.n_edges; j++)
        {
            Particle& pj = m_system->get_particle(vi.neigh[j]);
            Vertex& vj = mesh.get_vertices()[vi.neigh[j]];
            if (pj.in_tissue && (m_include_boundary || !vj.boundary))
            {
                if (m_has_part_params)
                {
                    K = m_particle_params[vj.type - 1].K;
                    gamma = m_particle_params[vj.type - 1].gamma;
                    lambda = m_particle_params[vj.type - 1].lambda;
                }
                double dA = 0.0;
                dA = (vj.area - pj.A0);
                double area_term = 0.5 * K * dA;
                double perim_term = gamma * vj.perim;
                double perim_term_edge_myosin = 0.0;
                double perim_term_perim_myosin = gamma * vj.perim
                        * (H_base + H_alpha * (pow(vj.perim / H_Ps, H_n)) / (H_KL + pow(vj.perim / H_Ps, H_n)));
                Vector3d area_vec(0.0, 0.0, 0.0);
                Vector3d perim_vec(0.0, 0.0, 0.0);
                Vector3d con_vec(0.0, 0.0, 0.0);
                Vector3d Nvec = Vector3d(pj.Nx, pj.Ny, pj.Nz);
                // My changes
                if (consider_edge_myosin_activity)
                {
                    for (int f = 0; f < vj.n_faces; f++)
                    {
                        int fid = vj.dual[f];
                        Face& f_nu = mesh.get_faces()[fid];
                        int fid_p = vj.dual[next(f, vj.n_faces)];
                        Face& f_nu_p = mesh.get_faces()[fid_p];
                        Vector3d& r_nu = f_nu.rc;
                        Vector3d& r_nu_p = f_nu_p.rc;
                        if (!(f_nu.is_hole || f_nu_p.is_hole))
                        {
                            double mij_p = H_base
                                    + H_alpha * pow((r_nu_p - r_nu).len() / H_Ls, H_n)
                                            / (H_KL + pow((r_nu_p - r_nu).len() / H_Ls, H_n));
                            perim_term_edge_myosin += gamma * sqrt(mij_p) * (r_nu_p - r_nu).len();
                        }
                    }
                }

                for (int f = 0; f < vj.n_faces; f++)
                {
                    int fid = vj.dual[f];
                    Face& f_nu = mesh.get_faces()[fid];
                    if (f_nu.has_vertex(i))
                    {
                        int fid_m = vj.dual[prev(f, vj.n_faces)];
                        int fid_p = vj.dual[next(f, vj.n_faces)];
                        Face& f_nu_m = mesh.get_faces()[fid_m];
                        Face& f_nu_p = mesh.get_faces()[fid_p];
                        Vector3d& r_nu_m = f_nu_m.rc;
                        Vector3d& r_nu = f_nu.rc;
                        Vector3d& r_nu_p = f_nu_p.rc;

                        Vector3d cross_prod_1(0.0, 0.0, 0.0);
                        if (!(f_nu.is_hole || f_nu_p.is_hole)) cross_prod_1 = cross(r_nu_p, Nvec)
                                * f_nu.get_jacobian(i);
                        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_1 = cross_prod_1
                                - cross(r_nu_m, Nvec) * f_nu.get_jacobian(i);
                        area_vec = area_vec + cross_prod_1;

                        // My changes
                        Vector3d cross_prod_2(0.0, 0.0, 0.0);
                        if (!(f_nu_m.is_hole || f_nu.is_hole))
                        {
                            if (consider_edge_myosin_activity)
                            {
                                double mij_m = H_base
                                        + H_alpha * pow((r_nu - r_nu_m).len() / H_Ls, H_n)
                                                / (H_KL + pow((r_nu - r_nu_m).len() / H_Ls, H_n));
                                cross_prod_2 = (sqrt(mij_m) * (r_nu - r_nu_m).unit()) * f_nu.get_jacobian(i);
                            }
                            else cross_prod_2 = ((r_nu - r_nu_m).unit()) * f_nu.get_jacobian(i);
                        }
                        if (!(f_nu_p.is_hole || f_nu.is_hole))
                        {
                            if (consider_edge_myosin_activity)
                            {
                                double mij_p = H_base
                                        + H_alpha * pow((r_nu_p - r_nu).len() / H_Ls, H_n)
                                                / (H_KL + pow((r_nu_p - r_nu).len() / H_Ls, H_n));
                                cross_prod_2 -= (sqrt(mij_p) * (r_nu_p - r_nu).unit()) * f_nu.get_jacobian(i);
                            }
                            else cross_prod_2 -= ((r_nu_p - r_nu).unit()) * f_nu.get_jacobian(i);
                        }

                        perim_vec = perim_vec + cross_prod_2;

                        // Adding force contributions from the dual edges associated with face f and face prev(f,vi.n_faces)
                        if (m_has_pair_params)
                        {
                            Vertex& vn = mesh.get_vertices()[vj.dual_neighbour_map[f]];
                            lambda = m_pair_params[vj.type - 1][vn.type - 1].lambda;
                        }
                        Vector3d cross_prod_3(0.0, 0.0, 0.0);
                        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_3 = lambda
                                * (((r_nu - r_nu_m).unit()) * f_nu.get_jacobian(i));
                        con_vec = con_vec + cross_prod_3;

                        // Adding force contributions from the dual edges associated with face f and face next(f,vi.n_faces)
                        if (m_has_pair_params)
                        {
                            Vertex& vn = mesh.get_vertices()[vj.dual_neighbour_map[next(f, vj.n_faces)]];
                            lambda = m_pair_params[vj.type - 1][vn.type - 1].lambda;
                        }
                        Vector3d cross_prod_4(0.0, 0.0, 0.0);
                        if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_4 = lambda
                                * (((r_nu_p - r_nu).unit()) * f_nu.get_jacobian(i));
                        con_vec = con_vec - cross_prod_4;

                    }
                }
                // area term
                double area_fact = -alpha * area_term;
                pi.fx += area_fact * area_vec.x;
                pi.fy += area_fact * area_vec.y;
                pi.fz += area_fact * area_vec.z;
                // perimeter term
                // My changes
                double perim_fact = -alpha
                        * (consider_edge_myosin_activity ? perim_term_edge_myosin :
                                (consider_perim_myosin_activity ? perim_term_perim_myosin : perim_term));

                pi.fx += perim_fact * perim_vec.x;
                pi.fy += perim_fact * perim_vec.y;
                pi.fz += perim_fact * perim_vec.z;
                // Add contractile term
                pi.fx -= alpha * con_vec.x;
                pi.fy -= alpha * con_vec.y;
                pi.fz -= alpha * con_vec.z;

                if (m_compute_stress)
                {
                    if (!vi.boundary && vi.area > 0)
                    {
                        double inv_area = 1.0 / vi.area;
                        pi.s_xx *= inv_area;
                        pi.s_xy *= inv_area;
                        pi.s_xz *= inv_area;
                        pi.s_yx *= inv_area;
                        pi.s_yy *= inv_area;
                        pi.s_yz *= inv_area;
                        pi.s_zx *= inv_area;
                        pi.s_zy *= inv_area;
                        pi.s_zz *= inv_area;
                    }
                }

            } // end of "if (pj.in_tissue && (m_include_boundary || !vj.boundary))" for neighboring cells
        } //end of checking neighbors

        // potential energy
        if (m_system->compute_per_particle_energy()) pi.add_pot_energy("my_vp", pot_eng);
    } // END of "for  (int i = 0; i < N; i++)";
} // end of compute(dt)
