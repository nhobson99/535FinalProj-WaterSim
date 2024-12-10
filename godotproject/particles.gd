extends Node

var SMOOTHING_RADIUS = 0.1


# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	# Initialize particles here if needed
	pass

# Smoothing function
func smoothing_function(d: float, radius: float) -> float:
	if d > radius:
		return 0.0
	var scale = 315.0 / (64.0 * PI * pow(radius, 9))
	var v = radius * radius - d * d
	return scale * v * v * v

# Density calculation
func density(v: Vector3, particles: Array) -> float:
	var d: float = 0.0
	for particle in particles:
		d += smoothing_function((v - particle.position).length(), SMOOTHING_RADIUS)
	return d

# Density gradient calculation
func density_gradient(v: Vector3, particles: Array) -> Vector3:
	var gradient = Vector3.ZERO
	for particle in particles:
		var diff = v - particle.position
		var distance = diff.length()
		if distance > SMOOTHING_RADIUS:
			continue
		var smoothing_value = smoothing_function(distance, SMOOTHING_RADIUS)
		gradient += diff.normalized() * smoothing_value
	return gradient

# Force calculation
func calc_force(p1: Vector3, p2: Vector3, d1: float, d2: float, grad1: Vector3, grad2: Vector3) -> Vector3:
	var diff = p1 - p2
	var distance = diff.length()
	
	if distance > SMOOTHING_RADIUS:
		return Vector3.ZERO
	
	var smoothing_value = smoothing_function(distance, SMOOTHING_RADIUS)
	var pressure_force = (d1 - d2) * smoothing_value * 0.001
	var gradient_force = (grad1 - grad2) * smoothing_value * 0.001
	
	return diff.normalized() * (pressure_force + gradient_force.length()) / d1


# Update function
func _process(delta: float) -> void:
	var children: Array[Node] = get_children()
	
	# Calculate densities and density gradients
	for particle in children:
		particle.set_meta("density", density(particle.position, children))
		particle.set_meta("density_gradient", density_gradient(particle.position, children))
	
	# Apply forces
	for particle: RigidBody3D in children:
		for other_particle in children:
			if particle == other_particle:
				continue
			var force = calc_force(
				particle.position, other_particle.position,
				particle.get_meta("density"), other_particle.get_meta("density"),
				particle.get_meta("density_gradient"), other_particle.get_meta("density_gradient")
			)
			particle.apply_force(force)
	
	# Enforce boundaries
	for particle: RigidBody3D in children:
		for axis in range(3):
			if particle.position[axis] > 1.0:
				particle.linear_velocity[axis] = -abs(particle.linear_velocity[axis]) * 0.1
				particle.position[axis] = 1.0
			elif particle.position[axis] < -1.0:
				particle.linear_velocity[axis] = abs(particle.linear_velocity[axis]) * 0.1
				particle.position[axis] = -1.0
		
