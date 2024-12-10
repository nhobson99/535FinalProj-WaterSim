extends Node3D

var prefab: PackedScene = null
var particles: Node = null
var framecount = 0

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	prefab = preload("res://sphere_prefab.tscn")
	particles = get_node(NodePath("../Particles"))
	for i in range(100):
		var new_obj: RigidBody3D = prefab.instantiate()
		particles.add_child(new_obj)
		new_obj.position = position + (Vector3(randf(), randf(), randf()) - Vector3(0.5, 0.5, 0.5))
		new_obj.linear_velocity = Vector3(0, -1, 0)
# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	if ((framecount % 1 == 0) and (particles.get_child_count() < 1024) and (1.0/delta >= 15.0)):
		var new_obj: RigidBody3D = prefab.instantiate()
		particles.add_child(new_obj)
		new_obj.position = position + Vector3(randf(), randf(), randf())*0.001
		new_obj.linear_velocity = Vector3(0, -1, 0)
	framecount += 1
