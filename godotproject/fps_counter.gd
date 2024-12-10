extends MeshInstance2D

var particles: Node = null

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	particles = get_node("/root/Node3D/ParticleSystem/Particles")

var counter = 0

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(_delta: float) -> void:
	counter += 1
	
	if counter % 20 == 0:
		mesh.text = "FPS: %.1f" % (Engine.get_frames_per_second())
		mesh.text += "\nParticles: %s" % (particles.get_child_count())
