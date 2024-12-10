shader_type canvas_item;
render_mode blend_mix;

// Uniforms for objects, camera, and settings
uniform vec3 camera_origin;
uniform vec3 camera_forward;
uniform vec3 camera_up;
uniform vec3 camera_right;
uniform float screen_distance;
uniform vec3[2048] sphere_pos;
uniform float[2048] sphere_velocities;
uniform float K: hint_range(0.001, 0.5);
uniform int max_steps: hint_range(1, 100);
uniform float threshold: hint_range(0.001, 0.1);
uniform int num_shapes: hint_range(0, 2048);
uniform float sphere_radius = 0.001f;

void vertex() {
	// Called for every vertex the material is visible on.
}

// Compute the signed distance to a sphere
float sdf_sphere(vec3 this_pos, vec3 other_pos, float radius) {
    return length(this_pos - other_pos) - radius;
}

float smoothMinPoly(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - h * (1.0 - h) * k;
}

struct best_hit {
	int index;
	float dist;
	float speed;
};

best_hit sdf(vec3 pos) {
	// Returns the index of the closest sphere
	best_hit output = best_hit(-1, 10000.0, 0.0);
	float mindist = 100.0;
    for (int i = 0; i < num_shapes; i++) {
        float dist_i = sdf_sphere(pos, sphere_pos[i], sphere_radius);
        mindist = smoothMinPoly(mindist, dist_i, K);
		if (mindist < output.dist) {
			output = best_hit(i, mindist, sphere_velocities[i]);
		}
    }
    return output;
}

// Perform raymarching
vec4 raymarch(vec3 ray_pos, vec3 ray_dir) {
	vec4 color = vec4(0);
    float total_distance = 0.0;
	int local_max_steps = max_steps;

    for (int step = 0; step < local_max_steps; step++) {
        vec3 point_on_ray = ray_pos + ray_dir * total_distance;
        best_hit hit = sdf(point_on_ray);
		if (hit.index == -1) continue;

        // If we're close enough to the surface, we hit it
        if (hit.dist < threshold) {
            color.a = 1.0;
			color.rgb += vec3(0.1, 0.1, 0.5) * 0.1;// + vec3(1.0, 1.0, 1.0) * hit.speed * 0.1;
			total_distance += sphere_radius * 2.0;
			local_max_steps = step + max_steps;
        }
		else {
        	total_distance += hit.dist;
		}
    }

    return color;
}

void fragment() {
	float aspect = SCREEN_PIXEL_SIZE.x / SCREEN_PIXEL_SIZE.y;
    vec2 screen_coord = (SCREEN_UV.xy - vec2(0.5, 0.5))*2.0;

	screen_coord.x /= aspect;

    // Calculate ray direction
    vec3 ray_dir = normalize(
		camera_forward * screen_distance
		+ camera_right * screen_coord.x
		- camera_up * screen_coord.y
	);

    // Perform raymarching from the camera origin along the ray direction
    vec4 raymarch_output = raymarch(camera_origin, ray_dir);

    // Set color based on distance (could add lighting here for shading)
    if (raymarch_output.a > 0.0) {
        COLOR = raymarch_output;
    } else {
        COLOR = vec4(0.0, 0.0, 0.0, 0.0);  // Background color
    }
}

//void light() {
	// Called for every pixel for every light affecting the material.
	// Uncomment to replace the default light processing function with this one.
//}
