extends ColorRect

var sphere_positions = []
var sphere_velocities = []
var particles: Node = null
var camera: Camera3D = null
var prefab: PackedScene = null

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	for i in range(2048):
		sphere_velocities.append(0.01)
	material.set_shader_parameter("sphere_pos", sphere_positions)
	material.set_shader_parameter("sphere_vel", sphere_velocities)
	particles = get_node(NodePath("/root/Node3D/ParticleSystem/Particles"))
	camera = get_node(NodePath("/root/Node3D/Camera3D"))
	prefab = preload("res://sphere_prefab.tscn")

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(_delta: float) -> void:
	sphere_positions = []
	sphere_velocities = []
	
	for particle: RigidBody3D in particles.get_children():
		sphere_positions.append(particle.position)
		sphere_velocities.append(particle.linear_velocity.length())
		
	#material.set_shader_parameter("sphere_radius", prefab.get_child(1).CSGMesh3D.Mesh.Radius)
	material.set_shader_parameter("sphere_pos", sphere_positions)
	material.set_shader_parameter("sphere_velocities", sphere_velocities)
	material.set_shader_parameter("num_shapes", sphere_positions.size())
	material.set_shader_parameter("camera_origin", camera.transform.origin)
	
	var camera_basis = camera.get_global_transform().basis
	material.set_shader_parameter("camera_forward", -camera_basis.z)
	material.set_shader_parameter("camera_up", camera_basis.y)
	material.set_shader_parameter("camera_right", camera_basis.x)
	
	
