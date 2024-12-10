extends RigidBody3D


# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	apply_force((Vector3(randf(), randf(), randf())-Vector3(0.5, 0.5, 0.5))*0.1)

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(_delta: float) -> void:
	#if abs(position.x) > 1.38 || abs(position.y) > 1.38:
		#self.queue_free()
	pass
