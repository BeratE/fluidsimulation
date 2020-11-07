#include "obj_reader.h"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include "tiny_obj_loader/tiny_obj_loader.h"


std::vector<learnSPH::TriMesh> learnSPH::readTriMeshesFromObj(std::string filename)
{
	std::vector<TriMesh> tri_meshes;

	// Initialise tinyobjloader objects and read the file
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string warn;
	std::string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str());

	// Input checking
	if (!warn.empty()) {
		std::cout << warn << std::endl;
	}

	if (!err.empty()) {
		std::cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	// Global information
	const int total_n_vertices = (int)attrib.vertices.size() / 3;

	// Write the geometric information into individual triangular meshes
	// Loop over meshes
	for (int shape_i = 0; shape_i < shapes.size(); shape_i++) {

		// Initialize individual triangular mesh
		TriMesh tri_mesh;
		tri_mesh.triangles.resize(shapes[shape_i].mesh.num_face_vertices.size());
		std::vector<bool> global_nodes_present(total_n_vertices, false);

		// Loop over triangles
		int index_offset = 0;
		for (int tri_i = 0; tri_i < shapes[shape_i].mesh.num_face_vertices.size(); tri_i++) {
			if (shapes[shape_i].mesh.num_face_vertices[tri_i] != 3) {
				std::cout << "learnSPH error: readTriMeshesFromObj can only read triangle meshes." << std::endl;
			}

			// Gather triangle global indices
			std::array<int, 3> triangle_global_indices;
			for (int vertex_i = 0; vertex_i < 3; vertex_i++) {
				tinyobj::index_t idx = shapes[shape_i].mesh.indices[(int)(3*tri_i + vertex_i)];
				const int global_vertex_index = idx.vertex_index;
				triangle_global_indices[vertex_i] = global_vertex_index;
				global_nodes_present[global_vertex_index] = true;
			}
			tri_mesh.triangles.push_back(triangle_global_indices);
		}

		// Reduce global indexes to local indexes
		std::vector<int> global_to_local_vertex_idx(total_n_vertices, -1);
		int local_vertices_count = 0;
		for (int global_vertex_i = 0; global_vertex_i < total_n_vertices; global_vertex_i++) {
			if (global_nodes_present[global_vertex_i]) {
				// Map global -> local
				global_to_local_vertex_idx[global_vertex_i] = local_vertices_count;
				local_vertices_count++;

				// Add vertex to the local mesh vertex vector
				tinyobj::real_t vx = attrib.vertices[(int)(3 * global_vertex_i + 0)];
				tinyobj::real_t vy = attrib.vertices[(int)(3 * global_vertex_i + 1)];
				tinyobj::real_t vz = attrib.vertices[(int)(3 * global_vertex_i + 2)];
				tri_mesh.vertices.push_back({vx, vy, vz});
			}
		}

		// Change triangle indices
		for (int tri_i = 0; tri_i < tri_mesh.triangles.size(); tri_i++) {
			for (int vertex_i = 0; vertex_i < 3; vertex_i++) {
				tri_mesh.triangles[tri_i][vertex_i] = global_to_local_vertex_idx[tri_mesh.triangles[tri_i][vertex_i]];
			}
		}

		tri_meshes.push_back(tri_mesh);
	}
	
	return tri_meshes;
}
