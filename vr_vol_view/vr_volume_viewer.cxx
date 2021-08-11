#include "vr_volume_viewer.h"
#include <cgv/math/ftransform.h>
#include <cgv/math/rigid_transform.h>

/// standard constructor
vr_volume_viewer::vr_volume_viewer() : plot("trigonometry", 2)
{
	translation_mode = false;
	rotation_mode = false;
	//deselecting_mode = false;
	scaling_mode = false;
	color_selected = false;
	trigger_button_down = false;
	show_slices_button_pressed = false;
	show_slices = false;
	translation_speed = 1.3;
	rotation_speed = 0.2f;
	trackball_size = 0.1f;
	index_of_currently_picked_sphere = -1;

	max_voxel_count_for_one_value = 0;
	max_voxel_count_for_one_value_in_second_histogram = 0;
	number_of_voxels_in_volume = 0;
	ray_length = 1;
	ray_color = rgb(0.5, 0.5, 0.5);
	range_of_neighboring_voxels_selected = 2;
	width_of_normal_distribution = 65;
	std_deviation = 8;
	texture_resolution = 256;
	normal_distribution_multiplier = 11;
	selected_color = vec3(1, 0, 0);
	position_relative_to_plot = vec3(0.6, -1, -0.25);
	color_picker_scaling = 0.5;

	for(auto& i : exact_slice_index)
	{
		i = 0;
	}
	vr_tf = dynamic_cast<vr_transfer_function*>(classif->get_transfer_functions().at(0));

	for (auto x = 0; x < texture_resolution; ++x)
	{
		hue_values_for_transfer_graph.push_back(1.01f);
		//hue_values_for_transfer_graph.push_back(P.at(x).x());
		hue_values_for_transfer_graph_example.push_back(1.01f);
	}
	hue_values_for_transfer_graph.push_back(0);
	hue_values_for_transfer_graph.push_back(1);
	hue_values_for_transfer_graph_example.push_back(0);
	hue_values_for_transfer_graph_example.push_back(1);
	
	plot.add_sub_plot("histogram");
	plot.add_sub_plot("selected_voxels_histogram");
	plot.add_sub_plot("transfer_function");
	plot.add_sub_plot("example_transfer_function");
	
	plot.set_extent(vecn(1.6f, 0.8f));
	plot.ref_sub_plot_samples(0).resize(texture_resolution);
	plot.ref_sub_plot_samples(1).resize(texture_resolution);
	plot.ref_sub_plot_samples(2).resize(texture_resolution);
	plot.ref_sub_plot_samples(3).resize(texture_resolution);
	for (auto x = 0; x < texture_resolution; ++x)
	{
		//plot.ref_sub_plot_samples(0)[x].set((float)x, 0);
		plot.ref_sub_plot_samples(1)[x].set(float(x) / float(texture_resolution), 0);
		plot.ref_sub_plot_samples(2)[x].set(float(x) / float(texture_resolution), float(x) / float(texture_resolution));
		plot.ref_sub_plot_samples(3)[x].set(float(x) / float(texture_resolution), float(x) / float(texture_resolution));
	}

	plot.set_sub_plot_attribute(2, 2, &hue_values_for_transfer_graph[0], hue_values_for_transfer_graph.size(), sizeof(float));
	plot.set_sub_plot_attribute(3, 2, &hue_values_for_transfer_graph_example[0], hue_values_for_transfer_graph_example.size(), sizeof(float));

	plot.ref_sub_plot2d_config(0).show_bars = true;
	plot.ref_sub_plot2d_config(0).bar_color.color = rgb(0.3f, 0.3f, 0.3f);
	plot.ref_sub_plot2d_config(0).bar_outline_color.color = rgb(0.6f, 0.6f, 0.6f);
	plot.ref_sub_plot2d_config(0).bar_percentual_width.size = 0.9f;
	plot.ref_sub_plot2d_config(0).bar_outline_width = 1.0f;
	plot.ref_sub_plot2d_config(0).show_lines = false;
	plot.ref_sub_plot2d_config(0).line_width.size = 2.0f;
	plot.ref_sub_plot2d_config(0).line_color.color = rgb(0.3f, 0.0f, 0.0f);
	plot.ref_sub_plot2d_config(0).show_points = false;
	plot.ref_sub_plot2d_config(1).show_bars = true;
	plot.ref_sub_plot2d_config(1).bar_color.color = rgb(0.0f, 0.0f, 1.0f);
	plot.ref_sub_plot2d_config(1).bar_percentual_width.size = 0.9f;
	plot.ref_sub_plot2d_config(1).bar_outline_width = 1.0f;
	plot.ref_sub_plot2d_config(1).show_lines = false;
	plot.ref_sub_plot2d_config(1).line_width.size = 2.0f;
	plot.ref_sub_plot2d_config(1).line_color.color = rgb(0.0f, 0.0f, 1.0f);
	plot.ref_sub_plot2d_config(1).show_points = false;
	//plot.ref_sub_plot2d_config(2).line_color.color = rgb(0.0f, 1.0f, 0.0f);
	plot.ref_sub_plot2d_config(3).show_plot = false;
	
	plot.color_mapping[0] = 2;
	plot.color_scale_index[0] = cgv::media::CS_HUE;
	plot.ref_sub_plot2d_config(2).set_color_indices(0);
	plot.ref_sub_plot2d_config(3).set_color_indices(0);

	plot.ensure_font_names();
	plot.set_label_font(16, cgv::media::font::FFA_BOLD_ITALIC, "Arial");

	// adjust domain, tick marks and extent in world space (of offline rendering process)
	plot.adjust_domain_to_data();
	plot.adjust_tick_marks();
	plot.adjust_extent_to_domain_aspect_ratio();
	plot.place_origin(vec3(-1, 1, 0));

	
	sphere_positions.emplace_back(vec3(0, 0, 0));
	sphere_positions.emplace_back(vec3(1, 0, 0));
	sphere_positions.emplace_back(vec3(2, 0, 0));

	sphere_colors.emplace_back(0, 0, 1);
	sphere_colors.emplace_back(1, 0, 0);
	sphere_colors.emplace_back(0, 1, 0);

	sphere_render_style.radius = 0.05f;

	set_hue_and_plot_values_for_example_graph();
}

float vr_volume_viewer::get_hsv_from_rgb(vec3 rgb)
{
	float max = std::max(std::max(rgb(0), rgb(1)), rgb(2));
	float min = std::min(std::min(rgb(0), rgb(1)), rgb(2));
	float delta = max - min;
	auto temp = 0.0;
	if(delta == 0)
	{
		return 0;
	}
	else if(max == rgb(0))
	{
		temp = rgb.y() - rgb.z();
		temp /= delta;
		temp = fmod(temp,6);
		temp /= 6;
	}
	else if (max == rgb(1))
	{
		temp = rgb.z() - rgb.x();
		temp /= delta;
		temp += 2;
		temp /= 6;
	}
	else if (max == rgb(2))
	{
		temp = rgb.x() - rgb.y();
		temp /= delta;
		temp += 4;
		temp /= 6;
	}
	if (temp < 0) return temp + 1;
	return temp;
}

void vr_volume_viewer::get_rgb_from_hsv(const float hue, const float saturation, vec3& rgb) const
{
	// C = s
	auto const x = /*saturation **/ (1 - abs(fmod(hue / 60, 2) - 1));
	//float m = 1 - saturation;
	if (hue >= 0 && hue < 60)
	{
		rgb = vec3(saturation, x, 0);
	}
	else if (hue >= 60 && hue < 120)
	{
		rgb = vec3(x, saturation, 0);
	}
	else if (hue >= 120 && hue < 180)
	{
		rgb = vec3(0, saturation, x);
	}
	else if (hue >= 180 && hue < 240)
	{
		rgb = vec3(0, x, saturation);
	}
	else if (hue >= 240 && hue < 300)
	{
		rgb = vec3(x, 0, saturation);
	}
	else if (hue >= 300 && hue <= 360)
	{
		rgb = vec3(saturation, 0, x);
	}
	//rgb += vec3(m, m, m);
	//rgb *= 255;
}

bool vr_volume_viewer::get_rotation_quat_from_controller_translation(quat& q, const vec3 controller_translation, vec3 view_direction) const
{
	// zero rotation
	if (controller_translation == vec3(0))
	{
		q.zeros();
		q.w() = 1;
		return false;
	}
	view_direction *= trackball_size;

	auto rotation_axis = cross(view_direction, controller_translation);
	if (rotation_axis == vec3(0))
	{
		return false;
	}
	// how much to rotate around that rotation_axis
	auto t = (controller_translation - view_direction).sqr_length() / (2 * trackball_size);

	// avoid problems with out-of-control values
	if (t > 1) t = 1;
	if (t < -1) t = -1;
	// how much to rotate about axis
	auto const phi = 2 * asin(t) * rotation_speed;

	// get quaternion from rotation axis
	rotation_axis.normalize();
	q.x() = rotation_axis.x();
	q.y() = rotation_axis.y();
	q.z() = rotation_axis.z();
	q.x() *= sin(phi / 2);
	q.y() *= sin(phi / 2);
	q.z() *= sin(phi / 2);
	q.w() = cos(phi / 2);
	return true;
}

void vr_volume_viewer::change_selected_voxels_histogram(vec3 voxel_to_select)
{
	auto instance = volume_instances(current_instance);
	
	auto value = instance.volume->get_voxel_component<int>(voxel_to_select(0), voxel_to_select(1), voxel_to_select(2), 0);

	// check if value is already in selected_voxels_histogram
	auto it = std::find_if(selected_voxels_histogram.begin(), selected_voxels_histogram.end(), [&value](const std::pair<int, int>& element) {return element.first == value; });
	if (it == selected_voxels_histogram.end())
	{
		/*if (deselecting_mode)
		{
			return;
		}*/
		selected_voxels_histogram.emplace_back(std::make_pair(value, 1/**number_of_voxels_in_volume / 100000*/));	
	}
	else
	{
		/*if (deselecting_mode)
		{
			it->second -= 1/**number_of_voxels_in_volume / 100000*//*;
			if (it->second <= 0)
			{
				selected_voxels_histogram.erase(it);
			}
			return;
		}*/
		it->second += 1/**number_of_voxels_in_volume / 100000*/;
	}
}

void vr_volume_viewer::select_voxels()
{
	const auto instance = volume_instances(current_instance);
	vec3 transformed_ray_point_position;
	if(show_slices)
	{
		auto d = intersect_ray_with_slice_planes();
		transformed_ray_point_position = ray_origin + ray_direction * d;
		transformed_ray_point_position = instance.vox_from_world(transformed_ray_point_position);
	}
	else
	{
		transformed_ray_point_position = instance.vox_from_world(ray_point_position);
	}

	if (instance.is_inside_vox(transformed_ray_point_position))
	{
		auto voxel_position = vec3(int(floor(transformed_ray_point_position(0))), int(floor(transformed_ray_point_position(1))), int(floor(transformed_ray_point_position(2))));

		//coloring voxels in volume
		if (color_selected)
		{
			const auto x = instance.volume->get_voxel_component<int>(voxel_position(0), voxel_position(1), voxel_position(2), 0);
			add_normal_distribution(x, 0.01);
			classif->transfer_parameters_updated();
			return;
		}
		
		// performance boost? if selected voxel is already selected, dont try to add all nearby voxels
		const auto it = std::find(selected_voxels.begin(), selected_voxels.end(), voxel_position);
		if (it != selected_voxels.end() && trigger_button_down /* ||
			it == selected_voxels.end() && deselecting_mode*/)
		{
			return;
		}
		// go through all neighbouring voxels
		for (auto x = int(voxel_position.x()) - range_of_neighboring_voxels_selected; x < int(voxel_position.x()) + range_of_neighboring_voxels_selected; ++x)
		{
			if (x < 0)	x = 0;
			if (x > instance.volume->get_dimensions()(0)) break;
			for (auto y = int(voxel_position.y()) - range_of_neighboring_voxels_selected; y < int(voxel_position.y()) + range_of_neighboring_voxels_selected; ++y)
			{
				if (y < 0)	y = 0;
				if (y > instance.volume->get_dimensions()(1)) break;
				for (auto z = int(voxel_position.z()) - range_of_neighboring_voxels_selected; z < int(voxel_position.z()) + range_of_neighboring_voxels_selected; ++z)
				{
					if (z < 0)	z = 0;
					if (z > instance.volume->get_dimensions()(2)) break;
					auto voxel_to_select = vec3(x, y, z);
					if (!instance.is_inside_vox(voxel_to_select))
					{
						continue;
					}

					// check if voxel is not selected
					const auto it1 = std::find(selected_voxels.begin(), selected_voxels.end(), voxel_to_select);
					if (it1 == selected_voxels.end())
					{
						selected_voxels.push_back(voxel_to_select);
						change_selected_voxels_histogram(voxel_to_select);
					} 
					/*else if (deselecting_mode)
					{
						selected_voxels.erase(it1);
						change_selected_voxels_histogram(voxel_to_select);	
					}*/
				}
			}
		}
	}
}

float vr_volume_viewer::intersect_ray_with_slice_planes()
{
	auto instance = volume_instances(current_instance);
	// x-slice
	auto x_plane_normal =  instance.orientation * vec3(1, 0, 0);
	auto d_min = intersect_line_with_plane(sphere_positions.at(1), ray_origin, x_plane_normal, ray_direction);
	
	// y-slice
	auto y_plane_normal = instance.orientation * vec3(0, 1, 0);
	auto d = intersect_line_with_plane(sphere_positions.at(2), ray_origin, y_plane_normal, ray_direction);
	if ((d_min > d && d > 0) || isnan(d_min)) d_min = d;
	
	// z-slice
	auto z_plane_normal = instance.orientation * vec3(0, 0, 1);
	d = intersect_line_with_plane(sphere_positions.at(0), ray_origin, z_plane_normal, ray_direction);
	if ((d_min > d && d > 0) || isnan(d_min)) d_min = d;
	return d_min;
}

float vr_volume_viewer::intersect_line_with_plane(vec3 point_on_plane,vec3 point_on_line,vec3 plane_normal, const vec3 line_direction) const
{
	const auto n = dot((point_on_plane - point_on_line), plane_normal);
	const auto d = dot(line_direction, plane_normal);
	if (d == 0) return NAN;
	return n / d;
}

float vr_volume_viewer::get_distance_from_line_to_point(vec3 point, vec3 point_on_line, vec3 line_direction) const
{
	const auto distance_point_to_point_on_line = point_on_line - point;
	auto dot = cgv::math::dot(distance_point_to_point_on_line, line_direction);
	auto d = distance_point_to_point_on_line.sqr_length() * line_direction.sqr_length() - dot*dot;
	d /= line_direction.sqr_length();
	return sqrt(d);
}


template <typename T>
cgv::type::uint32_type extract_uint32(const T* ptr) {
	return *ptr;
}

void vr_volume_viewer::on_volume_change(volume_instance& I, volume_instance_ext& J)
{
	// reset plots
	plot.ref_sub_plot_samples(0).clear();
	for (auto x = 0; x < texture_resolution; ++x)
	{
		plot.ref_sub_plot_samples(1)[x].set(float(x) / float(texture_resolution), 0);
	}
	
	if (current_instance != -1)
	{
		const auto volume = *I.volume;
		std::vector<cgv::type::uint32_type> histogram;
		const auto bytes_per_component = cgv::type::info::get_type_size(volume.get_component_type());
		if (histogram.empty()) {
			const auto nr_voxels = volume.get_nr_voxels();
			const auto nr_components = volume.get_nr_components();
			const auto bytes_per_voxel = nr_components * bytes_per_component;
			const cgv::type::uint8_type* ptr = volume.get_data_ptr<cgv::type::uint8_type>();

			histogram.resize(/*J.max_density - J.min_density + 1*/texture_resolution);
			std::fill(histogram.begin(), histogram.end(), 0);

			for (size_t j = 0; j < nr_voxels; ++j) {
				++histogram[extract_uint32(ptr) - min_density];
				ptr += bytes_per_voxel;
			}
		}
		plot.ref_sub_plot_samples(0).resize(texture_resolution);
		max_voxel_count_for_one_value = *std::max_element(std::begin(histogram), std::end(histogram));
		for (auto x = 0; x < texture_resolution; ++x)
		{
			plot.ref_sub_plot_samples(0)[x].set(float(x)/float(texture_resolution), float(histogram[x])/float(max_voxel_count_for_one_value));
		}
		plot.adjust_domain_to_data();
		plot.adjust_tick_marks();
		
		const auto min = min_value(I.volume->ref_extent()) * 2;
		I.volume->ref_extent().x() /= min;
		I.volume->ref_extent().y() /= min;
		I.volume->ref_extent().z() /= min;
		I.render_scale = 50.0f / I.volume->get_extent().length();
		number_of_voxels_in_volume = volume.get_nr_voxels();
		I.position = vec3(-2, 1, 0);
		for (auto i = 0; i < 3; ++i)
		{
			I.slice_index[i] = 0;
		}
	}
	volume_viewer::on_volume_change(I, J);
}

void vr_volume_viewer::change_color_value(float* current_value, const double target_value, const double factor) const
{
	const auto diff = abs(*current_value - target_value);
	if (diff < 0.01) *current_value = target_value;
	else if (*current_value < target_value) *current_value += diff / 20 * factor;
	else *current_value -= diff / 20 * factor;
}


void vr_volume_viewer::add_normal_distribution(const int x, const double add)
{
	for (auto current_x = x - width_of_normal_distribution / 2; current_x < x + width_of_normal_distribution/2; ++current_x)
	{
		// don't want to add values outside of the graph
		if (current_x >= 0 && current_x < texture_resolution)
		{
			// calculate new value by adding/subtracting normal distribution
			const auto t = float(x - current_x) / std_deviation;
			auto value_added = 1.0 / (sqrt(2 * M_PI) * std_deviation) * exp(-0.5f * t * t);
			value_added *= normal_distribution_multiplier * abs(add);
			auto new_value = plot.ref_sub_plot_samples(2)[current_x].y();
			if (add > 0 || color_selected)
			{
				new_value += value_added;
			}
			else
			{
				new_value -= value_added;
			}
			if (new_value < 0) new_value = 0;
			if (new_value > 1) new_value = 1;

			// if color selected update transfer graph and transfer function 
			if(color_selected)
			{
				// if it is the first color drawing in the transfer graph, set color values directly
				if (hue_values_for_transfer_graph.at(current_x) == 1.01f)
				{
					vr_tf->samples[current_x][0] = selected_color(0);
					vr_tf->samples[current_x][1] = selected_color(1);
					vr_tf->samples[current_x][2] = selected_color(2);
					
				}
				// update color values in samples in vr_transfer_function
				change_color_value(&vr_tf->samples[current_x][0], selected_color(0), abs(value_added * 100));
				change_color_value(&vr_tf->samples[current_x][1], selected_color(1), abs(value_added * 100));
				change_color_value(&vr_tf->samples[current_x][2], selected_color(2), abs(value_added * 100));
				// update hue value for transfer graph in plot 
				hue_values_for_transfer_graph.at(current_x) = get_hsv_from_rgb(vec3(vr_tf->samples[current_x][0], vr_tf->samples[current_x][1], vr_tf->samples[current_x][2]));
			}
			else
			{
				// update alpha value
				vr_tf->samples[current_x][3] = new_value;
				// update plot
				plot.ref_sub_plot_samples(2)[current_x].set(float(current_x) / float(texture_resolution), new_value);
			}		
		}
	}
}

void vr_volume_viewer::handle_left_controller_events(cgv::gui::vr_pose_event vrpe)
{
	auto& instance = volume_instances(current_instance);
	if (vrpe.get_trackable_index() == 0)
	{
		// update ray position
		vrpe.get_state().controller[0].put_ray(&ray_origin(0), &ray_direction(0));
		ray_direction.normalize();
		ray_point_position = ray_origin + ray_direction * ray_length;

		// draw in plot or move slices
		if (trigger_button_down)
		{
			// move slice index of currently picked sphere
			if (index_of_currently_picked_sphere != -1)
			{
				auto controller_to_volume_transformation = instance.orientation * (vrpe.get_position() - vrpe.get_last_position());
				auto max_slice_index = instance.volume->get_dimensions()(index_of_currently_picked_sphere);
				// first calculate the exact_slice_index(type double) and then assign it to the actual slice index(int)
				exact_slice_index[index_of_currently_picked_sphere] += controller_to_volume_transformation(index_of_currently_picked_sphere) / instance.get_box_vol().get_extent()(index_of_currently_picked_sphere) * max_slice_index;
				instance.slice_index[index_of_currently_picked_sphere] = exact_slice_index[index_of_currently_picked_sphere];
				// dont move slice index out of bounds
				if (instance.slice_index[index_of_currently_picked_sphere] > max_slice_index) instance.slice_index[index_of_currently_picked_sphere] = max_slice_index;
				if (instance.slice_index[index_of_currently_picked_sphere] < 0) instance.slice_index[index_of_currently_picked_sphere] = 0;

				// dont interact with anything else because controller ray could also intersect with other UI elements
				return;
			}

			// draw in plot to add a normal distribution on x location to change transfer function
			const auto plane_normal = cross(plot.get_corner(0) - plot.get_center(), plot.get_corner(2) - plot.get_center());
			const auto d = intersect_line_with_plane(plot.get_center(), ray_origin, plane_normal, ray_direction);
			if (!isnan(d))
			{
				const auto intersection_point = ray_origin + ray_direction * d;
				const auto selected_plot_position = (intersection_point - plot.get_origin()) * plot.get_orientation().get_matrix();
				const auto x = selected_plot_position.x() * float(texture_resolution) / plot.get_extent().x();
				const auto y = selected_plot_position.y() * float(texture_resolution) / plot.get_extent().y();
				if (x >= 0 && x < float(texture_resolution) && y >= 0 && y < float(texture_resolution))
				{
					auto const add = (vrpe.get_position() - vrpe.get_last_position()).y();
					add_normal_distribution(x, add);
					classif->transfer_parameters_updated();
				}
			}
			plot.adjust_domain_to_data();
			plot.adjust_tick_marks();
		}

		// draw in histogram
		if (trigger_button_down /* || deselecting_mode*/)
		{
			// select voxels in volume that then will be shown in histogram
			select_voxels();
			for (const auto pair : selected_voxels_histogram)
			{
				if (pair.second > max_voxel_count_for_one_value_in_second_histogram)
				{
					max_voxel_count_for_one_value_in_second_histogram = pair.second;
				}
			}
			
			for (auto& i : selected_voxels_histogram)
			{
				const auto x = i.first;
				if(plot.ref_sub_plot2d_config(0).show_plot)
				{
					plot.ref_sub_plot_samples(1)[x].set(float(x) / float(texture_resolution), float(i.second) / float(max_voxel_count_for_one_value));
				}
				else
				{
					plot.ref_sub_plot_samples(1)[x].set(float(x) / float(texture_resolution), float(i.second) / float(max_voxel_count_for_one_value_in_second_histogram));
				}
			}
		}
	}
}


void vr_volume_viewer::handle_pose_events(cgv::gui::event& e)
{
	const auto vrpe = dynamic_cast<cgv::gui::vr_pose_event&>(e);
	auto& instance = volume_instances(current_instance);

	// adjust plot orientation to hmd orientation
	const auto hmd_orientation = reinterpret_cast<const cgv::math::fmat<float, 3, 3>&>(vrpe.get_state().hmd.pose[0]);
	const auto hmd_translation = reinterpret_cast<const vec3&>(vrpe.get_state().hmd.pose[9]);
	plot.set_orientation(hmd_orientation);
	//plot.place_origin(hmd_translation +  hmd_orientation * vec3(-1, 0, -2));


	
	// get left controller transformations
	handle_left_controller_events(vrpe);

	//get right controller transformations
	if (vrpe.get_trackable_index() == 1)
	{
		if (translation_mode)
		{
			//I.position = vrpe.get_rotation_matrix() * (I.position - vrpe.get_last_position()) + vrpe.get_position(); // rotates volume with hmd and way smoother
			instance.position += (vrpe.get_position() - vrpe.get_last_position()) * translation_speed; // also works
			//instance.orientation *= vrpe.get_rotation_matrix();
		}
		if (!rotation_mode)
		{
			return;
		}

		// calculate controller translation to volume rotation
		const auto controller_translation = vrpe.get_position() - vrpe.get_last_position();
		if (controller_translation == vec3(0))
		{
			return;
		}
		// get view direction of user through controller orientation
		auto view_direction = vec3(0, 0, 1);
		view_direction = vrpe.get_orientation() * vec3(view_direction);
		view_direction.normalize();

		quat rotation_from_controller_translation;
		if (get_rotation_quat_from_controller_translation(rotation_from_controller_translation, controller_translation, view_direction))
		{
			//instance.orientation = rotation_from_controller_translation.get_matrix() * instance.orientation;
		}
		instance.orientation = vrpe.get_rotation_matrix() * instance.orientation;
	
		//instance.orientation = vrpe.get_orientation();
	}
}

bool vr_volume_viewer::handle_key_events(cgv::gui::event& e)
{
	auto& instance = volume_instances(current_instance);
	auto& vrke = static_cast<cgv::gui::vr_key_event&>(e);
	if (vrke.get_key() == 0) return false;
	switch (vrke.get_key()) {
	case vr::VR_MENU:
		if (vrke.get_controller_index() == 1)
		{
			show_slices_button_pressed = !show_slices_button_pressed;
			if(show_slices_button_pressed)
			{
				if(instance.show_volume == ShowMode::SM_HIDE)
				{
					instance.show_volume = ShowMode::SM_CLIP;
				}
				else
				{
					instance.show_volume = ShowMode::SM_HIDE;
				}
				show_slices = !show_slices;
			}
		}
		if(vrke.get_controller_index() == 0)
		{
			if(vrke.get_action() == cgv::gui::KA_PRESS)
			{
				plot.ref_sub_plot2d_config(0).show_plot = !plot.ref_sub_plot2d_config(0).show_plot;

				if(plot.ref_sub_plot2d_config(0).show_plot)
				{
					if (max_voxel_count_for_one_value_in_second_histogram == 0) break;
					const auto ratio = double(max_voxel_count_for_one_value_in_second_histogram) / max_voxel_count_for_one_value;
					for (auto& point : plot.ref_sub_plot_samples(1))
					{
						point.y() *= ratio;
					}
				}
				else
				{
					if (max_voxel_count_for_one_value_in_second_histogram == 0) break;
					const double ratio = double(max_voxel_count_for_one_value) / max_voxel_count_for_one_value_in_second_histogram;
					for (auto& point : plot.ref_sub_plot_samples(1))
					{
						point.y() *= ratio;
					}
				}
			}
		}
		return true;
	case vr::VR_GRIP:

		if (vrke.get_controller_index() == 1)
		{
			translation_mode = !translation_mode;
		}
		if (vrke.get_controller_index() == 0)
		{
			//deselecting_mode = !deselecting_mode;
			selected_voxels_histogram.clear();
			for (auto& point : plot.ref_sub_plot_samples(1))
			{
				point.y() = 0;
			}
		}
		return true;
	case vr::VR_DPAD_DOWN_LEFT:
		// get called when pressing touchpad in this location(down left) and when releasing it; 
		return true;
	case vr::VR_DPAD_DOWN:
		if (translation_mode && translation_speed > 0.1)
		{
			translation_speed -= 0.1f;
		}
		if (rotation_mode && rotation_speed > 0.01)
		{
			rotation_speed -= 0.01f;
		}
		if(vrke.get_controller_index() == 1 && !translation_mode && !rotation_mode)
		{
			instance.volume->ref_extent().x() *= 0.9;
			instance.volume->ref_extent().y() *= 0.9;
			instance.volume->ref_extent().z() *= 0.9;
			instance.render_scale = 50.0f / instance.volume->get_extent().length();
		}
		if (vrke.get_controller_index() == 0)
		{
			if (trigger_button_down)
			{
				normal_distribution_multiplier -= 0.1f * normal_distribution_multiplier;
			}
			else
			{
				ray_length -= 0.1;
			}
		}
		return true;
	case vr::VR_DPAD_DOWN_RIGHT:
		return true;
	case vr::VR_DPAD_LEFT:
		return true;
	case vr::VR_DPAD_RIGHT:
		return true;
	case vr::VR_DPAD_UP_LEFT:
		return true;
	case vr::VR_DPAD_UP:
		if (translation_mode)
		{
			translation_speed += 0.1;
		}
		if (rotation_mode)
		{
			rotation_speed += 0.1;
		}
		if (vrke.get_controller_index() == 1 && !translation_mode && !rotation_mode)
		{
			instance.volume->ref_extent().x() *= 1.1;
			instance.volume->ref_extent().y() *= 1.1;
			instance.volume->ref_extent().z() *= 1.1;
			instance.render_scale = 50.0f / instance.volume->get_extent().length();
		}
		if (vrke.get_controller_index() == 0)
		{
			if(trigger_button_down)
			{
				normal_distribution_multiplier += 0.1f * normal_distribution_multiplier;
			}
			else
			{
				ray_length += 0.1;
			}
		}
		return true;
	case vr::VR_DPAD_UP_RIGHT:
		return true;
	case vr::VR_INPUT0_TOUCH: // called when touching touchpad and when releasing
		return true;
	case vr::VR_INPUT0: // called when pressing in the middle of the touchpad and on release
		return true;
	case vr::VR_INPUT1_TOUCH: //trigger button: get called when pressing and when releasing
		return true;
	case vr::VR_INPUT1: //trigger button: gets called when holding; if holding first VR_INPUT1_TOUCH -> VR_INPUT1 -> release -> VR_INPUT1_TOUCH -> VR_INPUT1
		if (vrke.get_controller_index() == 1)
		{
			rotation_mode = !rotation_mode;
		}
		if (vrke.get_controller_index() == 0)
		{
			trigger_button_down = !trigger_button_down;
			/*if (trigger_button_down)
			{
				plot.get_domain_config_ptr()->axis_configs[1].set_log_scale(true);
			}
			else plot.get_domain_config_ptr()->axis_configs[1].set_log_scale(false);*/

			if (trigger_button_down)
			{
				// get slice through spheres
				for (auto i = 0; i < sphere_positions.size(); ++i)
				{
					if ((ray_origin - sphere_positions.at(i)).length() < 0.1)
					{
						if (i == 0) index_of_currently_picked_sphere = 2;
						else index_of_currently_picked_sphere = i - 1;
					}
				}

				// get color from color picker
				auto transformation_matrix = cgv::math::translate4<double>(plot.get_center() + plot.get_orientation().get_matrix() * position_relative_to_plot)
					* cgv::math::scale4<double>(dvec3(1, 1, 1) * color_picker_scaling)
					* plot.get_orientation().get_homogeneous_matrix();
				// get homogeneous representation
				cgv::math::fvec<double, 4> homo = vec4(color_picker.get_positions().at(0).x(), color_picker.get_positions().at(0).y(), color_picker.get_positions().at(0).z(), 1);
				auto color_picker_center_homo = transformation_matrix * homo;
				auto color_picker_center = vec3(color_picker_center_homo.x(), color_picker_center_homo.y(), color_picker_center_homo.z());
				//auto color_picker_center = color_picker.get_positions().at(0);
				//color_picker_center.x() *= -1;
				//color_picker_center.z() *= -1;
				auto color_picker_normal = plot.get_orientation().get_matrix() * vec3(0, 0, -1);
				color_picker_normal.to_vec().normalize();
				auto d = intersect_line_with_plane(color_picker_center, ray_origin, color_picker_normal, ray_direction);
				if (!isnan(d))
				{
					const auto intersection_point = ray_origin + d * ray_direction;
					const auto selected_color_picker_position = (intersection_point - color_picker_center) /** color_picker_orientation.get_matrix()*/;
					if ((selected_color_picker_position).length() <= color_picker_scaling)
					{
						auto hue = acos(dot(selected_color_picker_position, plot.get_orientation().get_matrix() * vec3(0, -1, 0)) / selected_color_picker_position.length());
						hue *= 180 / M_PI;
						if ((plot.get_orientation().inverse().get_matrix() * selected_color_picker_position).x() >= 0) hue = 360 - hue;
						//auto saturation = selected_color_picker_position.length();
						get_rgb_from_hsv(hue, 1/*saturation*/, selected_color);
						ray_color = rgb(selected_color.x(), selected_color.y(), selected_color.z());

						color_selected = true;
						return true;
					}
				}

				// check if controller ray intersects with slices
				vec3 transformed_ray_point_position;
				if (show_slices)
				{
					auto d = intersect_ray_with_slice_planes();
					//if (d < 0) return;
					transformed_ray_point_position = ray_origin + ray_direction * d;
					transformed_ray_point_position = instance.vox_from_world(transformed_ray_point_position);
				}
				else
				{
					transformed_ray_point_position = instance.vox_from_world(ray_point_position);
				}
				if (!instance.is_inside_vox(transformed_ray_point_position))
				{
					const auto plane_normal = cross(plot.get_corner(0) - plot.get_center(), plot.get_corner(2) - plot.get_center());
					d = intersect_line_with_plane(plot.get_center(), ray_origin, plane_normal, ray_direction);
					if (!isnan(d))
					{
						const auto intersection_point = ray_origin + ray_direction * d;
						const auto selected_plot_position = (intersection_point - plot.get_origin()) * plot.get_orientation().get_matrix();
						const auto x = selected_plot_position.x() * float(texture_resolution) / plot.get_extent().x();
						const auto y = selected_plot_position.y() * float(texture_resolution) / plot.get_extent().y();
						if (x < 0 || x > float(texture_resolution) || y < 0 || y > float(texture_resolution))
						{
							color_selected = false;
							ray_color = rgb(0.5, 0.5, 0.5);
						}
					}
					else
					{
						color_selected = false;
						ray_color = rgb(0.5, 0.5, 0.5);
					}
				}
			}
			else
			{
				index_of_currently_picked_sphere = -1;
			}
		}
		return true;
	default:
		std::cout << "Button pressed not found" << std::endl;
	}
	return false;
}

bool vr_volume_viewer::handle(cgv::gui::event& e)
{
	if (current_instance != -1)
	{
		switch (e.get_kind())
		{
			case cgv::gui::EID_POSE:
			{
				handle_pose_events(e);
			}
			case cgv::gui::EID_KEY:
			{
				return handle_key_events(e);
			}
			default:
			{
				std::cout << "event not found" << std::endl;
				return false;
			}
		}
	}
	return false;
}

///
void vr_volume_viewer::on_set(void* member_ptr)
{
	update_member(member_ptr);
	post_redraw();
	volume_viewer::on_set(member_ptr);
}
//
std::string vr_volume_viewer::get_type_name() const 
{
	return "vr_volume_viewer"; 
}

bool vr_volume_viewer::init(cgv::render::context& ctx)
{
	plot.set_view_ptr(find_view_as_node());
	plot.set_context(&ctx);


	if (color_picker.read(get_input_directory() + "/models/color_picker.obj")) 
	{
		const auto angle = 90.0 / 180.0 * M_PI;
		const auto rotation = quat(quat::X_AXIS, angle);
		color_picker.transform(rotation.get_matrix(), vec3(-1, 0, -1));
		
		color_picker_info.construct(ctx, color_picker);
		color_picker_info.bind(ctx, ctx.ref_surface_shader_program(true), true);
		auto& mats = color_picker_info.get_materials();
		if (!mats.empty()) 
		{
			//setup texture paths
			const auto di = mats[0]->add_image_file(get_input_directory() + "/images/ColorWheel1.jpg");

			//ensure the textures are loaded
			mats[0]->ensure_textures(ctx);

			//set texture indices
			mats[0]->set_diffuse_index(di);
		}

		/*ctx.mul_modelview_matrix(
			cgv::math::translate4<double>(vec3(0, 1.0, 0)) *
			cgv::math::scale4<double>(dvec3(1, 1, 1)) *
			cgv::math::rotate4<double>(90, vec3(1, 0, 0))
		);*/
	}

	return plot.init(ctx) && volume_viewer::init(ctx);
}

void vr_volume_viewer::init_frame(cgv::render::context& ctx)
{
	plot.init_frame(ctx);
	volume_viewer::init_frame(ctx);
}

void vr_volume_viewer::draw(cgv::render::context& ctx)
{
	// draw color_picker
	if (color_picker_info.is_constructed())
	{
		//push back a model view matrix to transform the rendering of this mesh
		ctx.push_modelview_matrix();

		// translate and scale
		double R = 0.5;
		//const auto hmd_orientation = reinterpret_cast<const cgv::math::fmat<float, 3, 3>&>(get_vr_kit(get_kit_ptr()).get_s);
		ctx.mul_modelview_matrix(cgv::math::translate4<double>(plot.get_center() + plot.get_orientation().get_matrix() * position_relative_to_plot)
			* cgv::math::scale4<double>(dvec3(1, 1, 1) * color_picker_scaling)
			* plot.get_orientation().get_homogeneous_matrix()
			//cgv::math::rotate4<double>(90, vec3(1, 0, 0))
		);

		//quat rotation = quat(cgv::render::render_types::quat::X_AXIS, 90);
		//color_picker.transform(rotation.get_matrix(), vec3(0, 5.0, 0));
		// 
		// actually draw the mesh
		color_picker_info.draw_all(ctx);

		// restore the previous transform
		ctx.pop_modelview_matrix();
	}

	plot.draw(ctx);

	// draw one ray from each controller
	const auto vr_view_ptr = dynamic_cast<vr_view_interactor*>(get_view_ptr());
	if (vr_view_ptr) {
		std::vector<vec3> P;
		std::vector<float> R;
		std::vector<rgb> C;
		const auto state_ptr = vr_view_ptr->get_current_vr_state();
		if (state_ptr)
		{
			for (auto ci = 0; ci < 2; ++ci)
			{
				if (state_ptr->controller[ci].status == vr::VRS_TRACKED)
				{
					vec3 ray_origin, ray_direction;
					state_ptr->controller[ci].put_ray(&ray_origin(0), &ray_direction(0));
					P.push_back(ray_origin);
					R.push_back(0.005f);
					P.push_back(ray_origin + ray_length * ray_direction);
					R.push_back(0.005f);
					rgb c(float(1 - ci), 0, float(ci));
					C.push_back(ray_color);
					C.push_back(ray_color);
				}
			}
		}
		if (!P.empty()) {
			auto& cr = cgv::render::ref_rounded_cone_renderer(ctx);
			//cr.set_render_style(cone_style);
			//cr.set_eye_position(vr_view_ptr->get_eye_of_kit());
			cr.set_position_array(ctx, P);
			cr.set_color_array(ctx, C);
			cr.set_radius_array(ctx, R);
			if (!cr.render(ctx, 0, P.size())) {
				auto& prog = ctx.ref_default_shader_program();
				const auto pi = prog.get_position_index();
				const auto ci = prog.get_color_index();
				cgv::render::attribute_array_binding::set_global_attribute_array(ctx, pi, P);
				cgv::render::attribute_array_binding::enable_global_array(ctx, pi);
				cgv::render::attribute_array_binding::set_global_attribute_array(ctx, ci, C);
				cgv::render::attribute_array_binding::enable_global_array(ctx, ci);
				glLineWidth(3);
				prog.enable(ctx);
				glDrawArrays(GL_LINES, 0, GLsizei(P.size()));
				prog.disable(ctx);
				cgv::render::attribute_array_binding::disable_global_array(ctx, pi);
				cgv::render::attribute_array_binding::disable_global_array(ctx, ci);
				glLineWidth(1);
			}
		}
	}

	// draw spheres for slice interaction
	if(show_slices)
	{
		const auto volume_position = volume_instances(current_instance).position;
		const auto volume_orientation = volume_instances(current_instance).orientation;
		auto volume_extent = volume_instances(current_instance).volume->get_box().get_extent();
		const auto slice_index_x = volume_instances(current_instance).slice_index[0];
		const auto slice_index_y = volume_instances(current_instance).slice_index[1];
		const auto slice_index_z = volume_instances(current_instance).slice_index[2];
		const auto max_slice_index_x = volume_instances(current_instance).volume->get_dimensions()(0);
		const auto max_slice_index_y = volume_instances(current_instance).volume->get_dimensions()(1);
		const auto max_slice_index_z = volume_instances(current_instance).volume->get_dimensions()(2);
		sphere_positions.at(0) = volume_position;
		sphere_positions.at(0) += volume_orientation * vec3(volume_extent.x() / 2, 0, (float(slice_index_z)/(float(max_slice_index_z)/2)-1)* volume_extent.z() / 2);
		sphere_positions.at(1) = volume_position;
		sphere_positions.at(1) += volume_orientation * vec3((float(slice_index_x) / (float(max_slice_index_x)/2) - 1) * volume_extent.x() / 2, volume_extent.y() / 2, 0);
		sphere_positions.at(2) = volume_position;
		sphere_positions.at(2) += volume_orientation * vec3(0, (float(slice_index_y) / (float(max_slice_index_y)/2) - 1) * volume_extent.y() / 2, volume_extent.z() / 2);
		auto& sphere_renderer = cgv::render::ref_sphere_renderer(ctx);
		sphere_renderer.set_position_array(ctx, sphere_positions);
		sphere_renderer.set_color_array(ctx, sphere_colors);
		sphere_renderer.set_render_style(sphere_render_style);
		sphere_renderer.render(ctx, 0, sphere_positions.size());
	}
	
	volume_viewer::draw(ctx);
}

/// draw textual information here
void vr_volume_viewer::finish_frame(cgv::render::context& ctx)
{
	volume_viewer::finish_frame(ctx);
}

///
void vr_volume_viewer::place_volume_instance_on_table(int instance_idx)
{
	// extract information to place instance
	auto& instance = volume_instances(instance_idx);
	const auto box = instance.get_box_vol();
	const auto scale = 1.0f / box.get_extent().length();
	auto height = 0.75f;
	if (get_scene_ptr())
		height = get_scene_ptr()->get_table_extent()[1];
	height += 0.51f * box.get_extent()[1] * scale;
	// place instance
	instance.position = vec3(0, height, 0);
	instance.volume->ref_extent() *= scale;
	instance.render_scale = 50.0f / instance.volume->get_extent().length();
	// ensure UI is up to date
	for (unsigned i = 0; i < 3; ++i) {
		on_set(&instance.position[i]);
		on_set(&instance.volume->ref_extent()[i]);
	}
	on_set(&instance.render_scale);
	// ensure render view is up to date
	post_redraw();
}

void vr_volume_viewer::set_hue_and_plot_values_for_example_graph()
{
	for(auto x = 0; x < texture_resolution; ++x)
	{
		auto n_x = float(x) / float(texture_resolution);
		auto f_value = ((n_x-0.6)*(n_x-0.3)*(n_x-0.7)*(n_x-0.5))/0.002 + 0.6;
		if (f_value > 1) f_value = 1;
		plot.ref_sub_plot_samples(3)[x].set(n_x, f_value);
		hue_values_for_transfer_graph_example.at(x) = f_value;
	}
}


void vr_volume_viewer::create_gui()
{
	add_member_control(this, "rotation_speed", rotation_speed, "value");
	add_decorator("vr volume viewer", "heading");
	connect_copy(add_button("place on table")->click, cgv::signal::rebind(this, &vr_volume_viewer::place_volume_instance_on_table, cgv::signal::_r(current_instance)));
	plot.create_gui(this, *this);
	volume_viewer::create_gui();
}

#ifndef NO_VR_VOLUME_VIEWER_INSTANCE
#include <cgv/base/register.h>
/// register a newly created cube with the name "cube1" as constructor argument
extern cgv::base::object_registration<vr_volume_viewer> vr_vol_viewer_reg("vr_volume_viewer");
#endif

#ifdef REGISTER_SHADER_FILES
#include <vol_view_shader_inc.h>
#endif

#ifdef CGV_FORCE_STATIC
extern cgv::base::registration_order_definition dro("stereo_view_interactor;vr_volume_viewer");
#endif


// Project an x, y pair onto a sphere of radius r OR a hyperbolic sheet if we are away from the center of the sphere.
/*float vr_volume_viewer::project_to_sphere(float r, float x, float y)
{
	float d, t, z;
	d = sqrt(x * x + y * y);
	// inside sphere
	if (d < r / std::sqrt(2))
	{
		z = sqrt(r * r - d * d);
	}
	// on hyperbola
	else
	{
		t = r / std::sqrt(2);
		z = t * t / d;
	}
	return z;
}*/