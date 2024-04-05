// From standard library
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

// eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// core
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier.hpp>
#include <minimesh/core/mohe/mesh_simplify.hpp>
#include <minimesh/core/mohe/mesh_deform.hpp>
#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/foldertools.hpp>
#include <minimesh/core/util/numbers.hpp>

// gui
#include <minimesh/viz/mesh_viewer.hpp>
#include <minimesh/viz/opengl_headers.hpp>


using namespace minimesh;

// ======================================================
// Global variables
// ======================================================
namespace globalvars
{
Mesh_viewer viewer;
mohe::Mesh_connectivity mesh;
mohe::Mesh_modifier modi(mesh);
mohe::Mesh_simplify simp(mesh);
mohe::Mesh_deform arap(mesh);
//
int glut_main_window_id;
//
GLUI * glui;
//
int num_entities_to_simplify;
bool initialized = false;
//
Eigen::Matrix3Xd displaced_vertex_positions;
}


// ======================================================
//              FREEGLUT CALL BACKS
// ======================================================
namespace freeglutcallback
{

void draw()
{
	globalvars::viewer.draw();
}


void window_reshaped(int w, int h)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.window_reshaped(w, h);

	if(should_redraw)
		glutPostRedisplay();
}


void keyboard_pressed(unsigned char c, int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.keyboard_pressed(c, x, y);
	// TODO IF THIS CHANGES WE SHOUD REINIT 
	if (c == 'f') {
		globalvars::arap.append_fixed(globalvars::arap.clickedVertex, globalvars::arap.clickedVertex);
		globalvars::arap.needs_reinit();
	}

	if (c == 'h') {
		globalvars::arap.append_handle(globalvars::arap.clickedVertex, globalvars::arap.clickedVertex);
		globalvars::arap.needs_reinit();
	}

	if(should_redraw)
		glutPostRedisplay();

}


void keyboard_arrows_pressed(int c, int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.keyboard_pressed(c, x, y);

	if(should_redraw)
		glutPostRedisplay();
}


void mouse_pushed(int button, int state, int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.mouse_pushed(button, state, x, y);

	//
	// NOTE: Sample of using Mesh_viewer for MESH DEFORMATION ASSINGMENT
	// Here is an example of how to use the selection feedback
	//
	{
		int clicked_on_vertex;
		bool did_user_click;
		globalvars::viewer.get_and_clear_vertex_selection(did_user_click, clicked_on_vertex);
		if (did_user_click) {
			globalvars::arap.clickedVertex = clicked_on_vertex;
			printf("User just clicked on vertex %d \n", clicked_on_vertex);
		}
	}

	if(should_redraw)
		glutPostRedisplay();
}


void mouse_moved(int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.mouse_moved(x, y);


	//
	// NOTE: Sample of using Mesh_viewer for MESH DEFORMATION ASSINGMENT
	// Here is an example of how to use the output from the viewer
	// If the user is displacing, I will displace the vertex being pulled
	//
	{
		bool has_pull_performed;
		Eigen::Vector3f pull_amount;
		int pulled_vert;
		globalvars::viewer.get_and_clear_vertex_displacement(has_pull_performed, pull_amount, pulled_vert);

		if(has_pull_performed)
		{
			force_assert(pulled_vert != Mesh_viewer::invalid_index);

			// globalvars::mesh.vertex_at(pulled_vert).data().xyz += pull_amount.cast<double>();

			if (globalvars::viewer.deform() == 1 && globalvars::arap.is_handle(pulled_vert)) {
				// printf("1");

				// deform the mesh
				globalvars::arap.deform(pulled_vert, pull_amount);

				// redraw
				{
					mohe::Mesh_connectivity::Defragmentation_maps defrag;
					globalvars::mesh.compute_defragmention_maps(defrag);
					globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
				}
			}

			// Must rerender now.
			should_redraw = true;
		}
	}

	if(should_redraw)
		glutPostRedisplay();
}


void subdivide_pressed(int)
{
	globalvars::modi.subdivision();

	// reload the mesh in the viewer
	mohe::Mesh_connectivity::Defragmentation_maps defrag;
	globalvars::mesh.compute_defragmention_maps(defrag);
	globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
	glutPostRedisplay();
}

void fixed_param_pressed(int)
{
	globalvars::modi.fixed_param();
    // mohe::Mesh_io(globalvars::mesh).write_obj("fixed_param.obj");

	// reload the mesh in the viewer
	mohe::Mesh_connectivity::Defragmentation_maps defrag;
	globalvars::mesh.compute_defragmention_maps(defrag);
	globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
	glutPostRedisplay();
}

void free_param_pressed(int)
{
	globalvars::modi.free_param();
    mohe::Mesh_io(globalvars::mesh).write_obj( "free_param.obj");

	// reload the mesh in the viewer
	mohe::Mesh_connectivity::Defragmentation_maps defrag;
	globalvars::mesh.compute_defragmention_maps(defrag);
	globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
	glutPostRedisplay();
}


void simplify_pressed(int)
{
	printf("Simplify button was pressed to remove %d entities \n", globalvars::num_entities_to_simplify);
	
	if (!globalvars::initialized) { 
		globalvars::simp.init();
		globalvars::initialized = true;
	} 
	globalvars::simp.simplify(globalvars::num_entities_to_simplify);

	// reload the mesh in the viewer
	mohe::Mesh_connectivity::Defragmentation_maps defrag;
	globalvars::mesh.compute_defragmention_maps(defrag);
	globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
	
	globalvars::simp.color_queue(globalvars::viewer.get_mesh_buffer(), defrag, globalvars::num_entities_to_simplify);
	glutPostRedisplay();
}


void show_spheres_pressed(int)
{

	// define anchor colors
	std::map<int, int> anchors = globalvars::arap.fixed_map();
	std::map<int, int> handles = globalvars::arap.handle_map();

	int a_size = static_cast<int>(anchors.size());
	int h_size = static_cast<int>(handles.size());

	Eigen::Matrix4Xf sphere_colors(4, a_size + h_size);
	Eigen::VectorXi sphere_indices(a_size + h_size);

	int i = 0;
	// update handle colors
	for (const auto &pair: anchors) {
		sphere_colors.col(i) << 1, 0, 0, 1;  // red for anchors
		sphere_indices(i) = pair.first;
		i++;
	}
	// update anchors colors
	for (const auto &pair: handles) {
		sphere_colors.col(i) << 0, 0, 1, 1;  // blue for handles
		sphere_indices(i) = pair.first;
		i++;
	}
	globalvars::viewer.get_mesh_buffer().set_colorful_spheres(sphere_indices, sphere_colors);

	glutPostRedisplay();

	
}

}

void clear_roi_pressed(int) {
	// clears the region of interest

	globalvars::arap.clear_constraints();
}

int main(int argc, char * argv[])
{
	// Remember current folder
	foldertools::pushd();

	// If no command line argument is specified, load a hardcoded mesh.
	// Useful when debugging with visual studio.
	// Change the hardcoded address to your needs.
	if(argc == 1)
	{
		// FOR MESHES W/O BOUNDARY
		foldertools::makeandsetdir("/Users/leofk/Documents/GitHub/modeling/mesh/");

		// A4
		mohe::Mesh_io(globalvars::mesh).read_auto("a4/woody-lo.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("a4/woody-hi.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("a4/bar.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("a4/cactus.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("a4/bumpy_plan.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("a4/cylinder.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("a4/hand.obj");
		

		// mohe::Mesh_io(globalvars::mesh).read_auto("cube.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("cow1.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("sphere1.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("camel.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("octopus.obj");

		// FOR MESHES W BOUNDARY
		// mohe::Mesh_io(globalvars::mesh).read_auto("pyramid.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("hexagon.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("cat.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("saddle.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("camel_head.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("mannequin.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("lion.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("balls.obj");
		// mohe::Mesh_io(globalvars::mesh).read_auto("hex_more.obj");
	}
	else // otherwise use the address specified in the command line
	{
		foldertools::makeandsetdir("/Users/leofk/Documents/524/modeling/mesh/");

		mohe::Mesh_io(globalvars::mesh).read_auto(argv[1]);
	}

	// Initialize GLUT window
	glutInit(&argc, argv);
	
	glutInitWindowSize(1920, 1080);
	glutInitDisplayMode(GLUT_STENCIL | GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE);
	globalvars::glut_main_window_id = glutCreateWindow("Mesh Viewer");

	// Initialize GLUI window for buttons and ...
	// globalvars::glui = GLUI_Master.create_glui("Controls");
	globalvars::glui = GLUI_Master.create_glui_subwindow(globalvars::glut_main_window_id, GLUI_SUBWINDOW_LEFT);
	globalvars::glui->set_main_gfx_window(globalvars::glut_main_window_id);

	// Register callbacks
	glutDisplayFunc(freeglutcallback::draw);
	GLUI_Master.set_glutReshapeFunc(freeglutcallback::window_reshaped);
	GLUI_Master.set_glutKeyboardFunc(freeglutcallback::keyboard_pressed);
	GLUI_Master.set_glutSpecialFunc(freeglutcallback::keyboard_arrows_pressed);
	GLUI_Master.set_glutMouseFunc(freeglutcallback::mouse_pushed);
	glutMotionFunc(freeglutcallback::mouse_moved);
	GLUI_Master.set_glutIdleFunc(NULL);

	// Initialize the viewer (it needs the bounding box of the mesh)
	Eigen::AlignedBox3f bbox;
	for(int v = 0; v < globalvars::mesh.n_total_vertices(); ++v)
	{
		mohe::Mesh_connectivity::Vertex_iterator vertex = globalvars::mesh.vertex_at(v);
		if(vertex.is_active())
		{
			bbox.extend(vertex.xyz().cast<float>());
		}
	}
	globalvars::viewer.initialize(bbox);
	
	// Load the mesh in the viewer
	{
		mohe::Mesh_connectivity::Defragmentation_maps defrag;
		globalvars::mesh.compute_defragmention_maps(defrag);
		globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
	}

	//
	// Add radio buttons to choose mesh
	//
	// GLUI_Panel * panel_mesh = globalvars::glui->add_panel("Choose Mesh");
	// GLUI_RadioGroup * radio_group_mesh = globalvars::glui->add_radiogroup_to_panel(panel_mesh, &globalvars::viewer.get_mesh());
	// for(int i = 0; i < Mesh_viewer::MESH_INVALID; ++i)
	// {
	// 	if(i == Mesh_viewer::MESH_CUBE)
	// 		globalvars::glui->add_radiobutton_to_group(radio_group_mesh, "Cube");
	// 	if(i == Mesh_viewer::MESH_COW)
	// 		globalvars::glui->add_radiobutton_to_group(radio_group_mesh, "Cow");
	// 	if(i == Mesh_viewer::MESH_PYRAMID)
	// 		globalvars::glui->add_radiobutton_to_group(radio_group_mesh, "Pyramid (w/ boundary)");
	// }

	//
	// Add radio buttons to see which mesh components to view
	// Please view GLUI's user manual to learn more.
	//

	GLUI_Panel * panel_view = globalvars::glui->add_panel("View mesh components");
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show vertices", &globalvars::viewer.get_draw_vertices());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show edges", &globalvars::viewer.get_draw_edges());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show faces", &globalvars::viewer.get_draw_faces());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show axis", &globalvars::viewer.get_draw_axis());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show lighting", &globalvars::viewer.get_has_lighting());

	//
	// Add radio buttons to determine mouse left click functionality
	//
	GLUI_Panel * panel_mouse_func = globalvars::glui->add_panel("Mouse functionality");
	GLUI_RadioGroup * radio_group_mouse_func =   globalvars::glui->add_radiogroup_to_panel(panel_mouse_func, &globalvars::viewer.get_mouse_function());
	for(int i = 0; i < Mesh_viewer::MOUSE_INVALID; ++i)
	{
		if(i == Mesh_viewer::MOUSE_VIEW)
			globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Pan and zoom");
		if(i == Mesh_viewer::MOUSE_SELECT)
			globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Select vertex");
		if(i == Mesh_viewer::MOUSE_MOVE_VERTEX)
			globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Move vertex");
	}

	//
	// Add subdivide button
	//
	GLUI_Button* button_subdivide =  globalvars::glui->add_button("Subdivide Loop", -1, freeglutcallback::subdivide_pressed);
	button_subdivide->set_w(200);

	GLUI_Button* button_fixed_param =  globalvars::glui->add_button("Harmonic Parameterization", -1, freeglutcallback::fixed_param_pressed);
	button_fixed_param->set_w(200);

	GLUI_Button* button_free_param =  globalvars::glui->add_button("LSCM Parameterization", -1, freeglutcallback::free_param_pressed);
	button_free_param->set_w(200);

	//
	// Add simplify button and a spinner to read how many entities to remove
	//
	globalvars::num_entities_to_simplify = 0;
	GLUI_Spinner* spinner_simplify = globalvars::glui->add_spinner("# of entities to simplify", GLUI_SPINNER_INT, &globalvars::num_entities_to_simplify);
	spinner_simplify->set_alignment(GLUI_ALIGN_CENTER);
	spinner_simplify->set_w(300);
	
	GLUI_Button* button_simplify = globalvars::glui->add_button("Simplify", -1, freeglutcallback::simplify_pressed);
	button_simplify->set_w(200);

	//
	// Add show spheres button to demo how to draw spheres on top of the vertices
	//
	globalvars::glui->add_button("Demo Showing Spheres", -1, freeglutcallback::show_spheres_pressed);

    //
    // Add ARAP Controls
    //

    // arap control options
    GLUI_Panel *arap_panel_options = globalvars::glui->add_panel("ARAP Deformation");
    globalvars::glui->add_checkbox_to_panel(arap_panel_options, "Deform", &globalvars::viewer.deform());

    // clear roi of interest TODO
    // GLUI_Button *clear_roi_button = globalvars::glui->add_button("Clear ROI", -1,
    //                                                              freeglutcallback::clear_roi_pressed);
    // clear_roi_button->set_w(200);


	//
	// Save the initial vertex positions
	//
	globalvars::displaced_vertex_positions.resize(3,globalvars::mesh.n_active_vertices() );
	for (int i = 0 ; i < globalvars::mesh.n_active_vertices() ; ++i)
	{
		globalvars::displaced_vertex_positions.col(i) = globalvars::mesh.vertex_at(i).xyz();
	}

	// Sync all glui variables
	globalvars::glui->sync_live();

	// Start main loop
	glutPostRedisplay(); // Draw everything again just for caution.
	glutMainLoop();

	// revert back to initial folder
	foldertools::popd();

	return 0;
}
