extends Camera3D

var origin: Vector3 = Vector3(0.0, -1.0, 0.0)
var phi: float = 0.0
var tau: float = 0.0
var radius: float = 1.0
var target_radius: float = 1.0
var smoothing_factor: float = 0.75

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	pass # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(_delta: float) -> void:
	radius = (radius * smoothing_factor) + (target_radius * (1.0 - smoothing_factor))
	position = origin + Vector3(
		radius * cos(tau) * cos(phi),
		radius * sin(tau),
		radius * cos(tau) * sin(phi)
	)
	look_at(origin)

var clicking: bool = false
var old_pos: Vector2 = Vector2(0.0, 0.0)

func _input(event):
	if event is InputEventMouseButton:
		if (event.button_index == MOUSE_BUTTON_LEFT
				or event.button_index == MOUSE_BUTTON_MIDDLE
				or event.button_index == MOUSE_BUTTON_RIGHT):
			clicking = event.pressed
			old_pos = event.position
		elif event.button_index == MOUSE_BUTTON_WHEEL_UP:
			target_radius -= 0.1
		elif event.button_index == MOUSE_BUTTON_WHEEL_DOWN:
			target_radius += 0.1
	elif event is InputEventMouseMotion and clicking:
		var diff = (event.position - old_pos)*0.01
		old_pos = event.position
		phi += diff.x
		tau += diff.y
	target_radius = min(3.0, max(0.5, target_radius))
	tau = min(PI*0.45, max(PI*0.05, tau))
	position = origin + Vector3(
		radius * cos(tau) * cos(phi),
		radius * sin(tau),
		radius * cos(tau) * sin(phi)
	)
	look_at(origin)
