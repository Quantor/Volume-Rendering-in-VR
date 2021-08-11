#pragma once

#include <cgv/defines/quote.h>
#include "../vol_view/volume_viewer.h"
#include "vol_viz/vr_transfer_function.h"
#include <plugins/vr_lab/vr_tool.h>
#include <cg_vr/vr_events.h>
#include "lib_begin.h"
#include "../libs/vol_viz/graph_data.h"


//template <class D, typename P = unsigned char*>
class CGV_API vr_volume_viewer :
	public volume_viewer,
	public vr::vr_tool
{
private:
	//true if button is pressed to translate volume; if true map controller translation to volume translation
	bool translation_mode;
	//true if button is pressed to rotate volume; if true map controller translation to volume rotation
	bool rotation_mode;
	//defines how much controller translation is transferred to volume translation in translation_mode
	float translation_speed;
	//defines how much controller translation is transferred to volume rotation in rotation_mode
	float rotation_speed;
	bool scaling_mode;
	bool color_selected;
	bool trigger_button_down;
	//bool deselecting_mode;
	bool show_slices;
	bool show_slices_button_pressed;

	float trackball_size;
	
	vec3 ray_point_position;
	vec3 ray_origin;
	vec3 ray_direction;
	
	float ray_length;
	rgb ray_color;
	int range_of_neighboring_voxels_selected;
	int width_of_normal_distribution;
	float std_deviation;
	int texture_resolution;
	int max_voxel_count_for_one_value;
	int max_voxel_count_for_one_value_in_second_histogram;
	int number_of_voxels_in_volume;
	float normal_distribution_multiplier;
	// vector to not select voxels twice
	std::vector<vec3> selected_voxels;
	// vector of pairs(first: voxel value, second: number of selected voxels with this value)
	std::vector<std::pair<int, int>> selected_voxels_histogram;
	std::vector<float> hue_values_for_transfer_graph;
	std::vector<float> hue_values_for_transfer_graph_example;
	
	vr_transfer_function* vr_tf;

	cgv::render::mesh_render_info color_picker_info;
	cgv::media::mesh::simple_mesh<> color_picker;
	vec3 selected_color;
	vec3 position_relative_to_plot;
	double color_picker_scaling;

	// index 0 -> z-slice, index 1 -> x-slice, index 2 -> y-slice
	std::vector<vec3> sphere_positions;
	std::vector<rgb> sphere_colors;
	cgv::render::sphere_render_style sphere_render_style;
	int index_of_currently_picked_sphere;
	double exact_slice_index[3]{};
	cgv::plot::plot2d plot;

protected:
	static std::string get_input_directory() { return QUOTE_SYMBOL_VALUE(INPUT_DIR); };

	float get_hsv_from_rgb(vec3 rgb);
	void get_rgb_from_hsv(float hue, float saturation, vec3& rgb) const;
	bool get_rotation_quat_from_controller_translation(quat& q, vec3 controller_translation, vec3 view_direction) const;
	// add/delete or increase/decrease voxel in selected_voxels_histogram
	void change_selected_voxels_histogram(vec3 voxel);

	// select range_of_neighboring_voxels_selected voxels around ray_point_position
	void select_voxels();

	float get_distance_from_line_to_point(vec3 point, vec3 point_on_line1, vec3 line_direction) const;
	float intersect_ray_with_slice_planes();
	float intersect_line_with_plane(vec3 point_on_plane, vec3 point_on_line, vec3 plane_normal, vec3 line_direction) const;
	void on_volume_change(volume_instance& I, volume_instance_ext& J) override;
	void change_color_value(float* current_value, double target_value, double factor) const;
	void add_normal_distribution(int x, double add);
	void handle_left_controller_events(cgv::gui::vr_pose_event vrpe);
	void handle_pose_events(cgv::gui::event& e);
	bool handle_key_events(cgv::gui::event& e);
	void place_volume_instance_on_table(int instance_idx);

	void set_hue_and_plot_values_for_example_graph();

public:
	/// standard constructor
	vr_volume_viewer();
	//float project_to_sphere(float r, float x, float y);

	bool handle(cgv::gui::event& e) override;
	void init_frame(cgv::render::context& ctx) override;
	bool init(cgv::render::context& ctx) override;

	///
	void on_set(void* member_ptr) override;
	//
	std::string get_type_name() const override;
	//
	void draw(cgv::render::context& ctx) override;
	/// draw textual information here
	void finish_frame(cgv::render::context&) override;
	void create_gui() override;
};

#include <cgv/config/lib_end.h>
