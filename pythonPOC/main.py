import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import multiprocessing
import time

plt.switch_backend("TkAgg")

pool = None

# Simulation parameters
NUM_PARTICLES = 100
DOMAIN_SIZE = 10.0
NUM_DOMAINS = 10
SMOOTHING_LENGTH = DOMAIN_SIZE / NUM_DOMAINS
TIME_STEP = 1.0/60
NUM_STEPS = 5000
PARTICLE_MASS = 1.0
STIFFNESS = 100.0  # Pressure coefficient
REST_DENSITY = 2.0  # Rest density for fluid
VISCOSITY = 0.05  # Viscosity coefficient
GRAVITY = 9.8

pressures = None

def spatial_hash(p: np.ndarray):
    global NUM_PARTICLES, DOMAIN_SIZE, SMOOTHING_LENGTH, TIME_STEP, \
        PARTICLE_MASS, STIFFNESS, REST_DENSITY, VISCOSITY, GRAVITY
    p0 = int(p[0] / SMOOTHING_LENGTH)
    p1 = int(p[1] / SMOOTHING_LENGTH)
    return p0*NUM_DOMAINS + p1

def hash_and_sort(particles):
    hashed_particles = [(particles[i], spatial_hash(particles[i])) for i in range(len(particles))]
    sorted_particles = sorted(hashed_particles, key=lambda x: x[1])

    particle_list = []
    indices = {}
    hashes = []
    for i in range(len(sorted_particles)):
        if (i == 0 or sorted_particles[i][1] != sorted_particles[i-1][1]):
            indices[sorted_particles[i][1]] = i
        particle_list.append(sorted_particles[i][0])
        hashes.append(sorted_particles[i][1])
        
    return np.asarray(particle_list), indices, hashes

def smoothing_kernel(r, h):
    """Cubic spline smoothing kernel."""
    q = r / h
    if q <= 1.0:
        return (1 / (np.pi * h**3)) * (1 - 1.5 * q**2 + 0.75 * q**3)
    elif q <= 2.0:
        return (1 / (np.pi * h**3)) * 0.25 * (2 - q)**3
    else:
        return 0.0

def smoothing_kernel_gradient(r_vec, r, h):
    """Gradient of the cubic spline smoothing kernel."""
    q = r / h
    r += (np.random.rand())*0.00001
    if q <= 1.0:
        return (-3 / (np.pi * h**4)) * (1 - q) * r_vec / r
    elif q <= 2.0:
        return (-3 / (4 * np.pi * h**4)) * (2 - q)**2 * r_vec / r
    else:
        return np.zeros_like(r_vec)

def surrounding_hashes(p: np.ndarray):
    global NUM_PARTICLES, DOMAIN_SIZE, SMOOTHING_LENGTH, TIME_STEP, \
        PARTICLE_MASS, STIFFNESS, REST_DENSITY, VISCOSITY, GRAVITY
    rv = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            rv.append(spatial_hash(p + np.array([i * SMOOTHING_LENGTH, j * SMOOTHING_LENGTH])))
    return rv

def compute_density(particles: np.ndarray, indices: dict, hashes: list, h):
    global NUM_PARTICLES, DOMAIN_SIZE, SMOOTHING_LENGTH, TIME_STEP, \
        PARTICLE_MASS, STIFFNESS, REST_DENSITY, VISCOSITY, GRAVITY
    """Compute densities for all particles."""
    densities = np.zeros(len(particles))
    for i, pi in enumerate(particles):
        hashes_to_search = surrounding_hashes(pi)
        for starting_hash in hashes_to_search:
            index = indices.get(starting_hash)
            if index == None:
                continue
            if index >= len(hashes):
                continue
            if (hashes[index] != starting_hash):
                continue
            while index < len(hashes) and hashes[index] == starting_hash:
                pj = particles[index]
                r = np.linalg.norm(pi - pj)
                densities[i] += PARTICLE_MASS * smoothing_kernel(r, h)
                index += 1
    return densities

def _compute_forces(args):
    """Compute forces for a single particle."""
    i, pi, particles, indices, hashes, densities, pressures, velocities, NUM_PARTICLES, DOMAIN_SIZE, SMOOTHING_LENGTH, TIME_STEP, \
        PARTICLE_MASS, STIFFNESS, REST_DENSITY, VISCOSITY, GRAVITY = args
    h = SMOOTHING_LENGTH
    force = np.zeros(2)
    hashes_to_search = surrounding_hashes(pi)
    for starting_hash in hashes_to_search:
            index = indices.get(starting_hash)
            if index == None:
                continue
            if index >= len(hashes):
                continue
            if (hashes[index] != starting_hash):
                continue
            while index < len(hashes) and hashes[index] == starting_hash:
                j, pj = index, particles[index]
                if i != j:
                    r_vec = pi - pj
                    r = np.linalg.norm(r_vec)
                    if r < 2 * h:
                        grad_w = smoothing_kernel_gradient(r_vec, r, h)

                        # Pressure force
                        pressure_term = -PARTICLE_MASS * (
                            pressures[i] / (densities[i]**2) + pressures[j] / (densities[j]**2))
                        force += pressure_term * grad_w

                        # Viscosity force
                        velocity_diff = velocities[j] - velocities[i]
                        viscosity_term = VISCOSITY * PARTICLE_MASS * velocity_diff / densities[j]
                        force += viscosity_term * smoothing_kernel(r, h)
                index += 1

    # Add gravitational force
    force += np.array([0, -GRAVITY]) * densities[i] + (np.random.rand(2) - 0.5) * 0.0001

    return force


def compute_forces(particles, indices, hashes, densities, velocities):
    """Compute pressure and viscosity forces for all particles."""
    global pressures
    global NUM_PARTICLES, DOMAIN_SIZE, SMOOTHING_LENGTH, TIME_STEP, \
        PARTICLE_MASS, STIFFNESS, REST_DENSITY, VISCOSITY, GRAVITY
    pressures = STIFFNESS * (densities - REST_DENSITY)
    # Prepare arguments for parallel processing
    args = [
        (i, pi, particles, indices, hashes, densities, pressures, velocities,
        NUM_PARTICLES, DOMAIN_SIZE, SMOOTHING_LENGTH, TIME_STEP,
        PARTICLE_MASS, STIFFNESS, REST_DENSITY, VISCOSITY, GRAVITY)
        for i, pi in enumerate(particles)
    ]
    
    forces = pool.map(_compute_forces, args)
    return np.array(forces)


def enforce_boundaries(particles, velocities, domain_size):
    """Ensure particles stay within the domain by reflecting them off boundaries."""
    for i, p in enumerate(particles):
        for d in range(len(p)):
            if p[d] < 0:
                particles[i][d] = 0
                velocities[i][d] = abs(velocities[i][d])
            elif p[d] > domain_size:
                particles[i][d] = domain_size
                velocities[i][d] = -abs(velocities[i][d])

def update_sliders(val):
    """Update simulation parameters based on slider values."""
    global NUM_PARTICLES, DOMAIN_SIZE, SMOOTHING_LENGTH, TIME_STEP, \
        PARTICLE_MASS, STIFFNESS, REST_DENSITY, VISCOSITY, GRAVITY

    NUM_PARTICLES = int(sliders["NUM_PARTICLES"].val)
    DOMAIN_SIZE = sliders["DOMAIN_SIZE"].val
    SMOOTHING_LENGTH = sliders["SMOOTHING_LENGTH"].val
    TIME_STEP = sliders["TIME_STEP"].val
    PARTICLE_MASS = sliders["PARTICLE_MASS"].val
    STIFFNESS = sliders["STIFFNESS"].val
    REST_DENSITY = sliders["REST_DENSITY"].val
    VISCOSITY = sliders["VISCOSITY"].val
    GRAVITY = sliders["GRAVITY"].val

def main():
    global velocities, densities, sliders, TIME_STEP
    # Initialize particle positions and velocities
    particles = np.random.rand(NUM_PARTICLES, 2) * DOMAIN_SIZE
    velocities = np.zeros((NUM_PARTICLES, 2))

    # Create a figure with subplots for sliders and visualization
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.set_box_aspect(1)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.subplots_adjust(left=0.25, right=0.75, bottom=0.4)
    ax.set_xlim(0, DOMAIN_SIZE)
    ax.set_ylim(0, DOMAIN_SIZE)

    scatter = ax.scatter([], [], c=[], cmap='viridis', s=200, vmin=0.0, vmax=REST_DENSITY * 2, alpha=0.5)
    plt.colorbar(scatter, label='Density')

    # Create slider axes and sliders
    slider_axes = [
        plt.axes([0.25, 0.35 - i * 0.04, 0.65, 0.03]) for i in range(9)
    ]
    sliders = {
        "NUM_PARTICLES": Slider(slider_axes[0], "Particles", 10, 1000, valinit=NUM_PARTICLES, valstep=10),
        "DOMAIN_SIZE": Slider(slider_axes[1], "Domain", 1.0, 50.0, valinit=DOMAIN_SIZE),
        "SMOOTHING_LENGTH": Slider(slider_axes[2], "Smooth", 0.1, 5.0, valinit=SMOOTHING_LENGTH),
        "TIME_STEP": Slider(slider_axes[3], "Time Step", 0.01, 1.0, valinit=TIME_STEP),
        "PARTICLE_MASS": Slider(slider_axes[4], "Mass", 0.1, 10.0, valinit=PARTICLE_MASS),
        "STIFFNESS": Slider(slider_axes[5], "Stiffness", 1.0, 200.0, valinit=STIFFNESS),
        "REST_DENSITY": Slider(slider_axes[6], "Rest Density", 0.1, 10.0, valinit=REST_DENSITY),
        "VISCOSITY": Slider(slider_axes[7], "Viscosity", 0.0, 1.0, valinit=VISCOSITY),
        "GRAVITY": Slider(slider_axes[8], "Gravity", 0.0, 20.0, valinit=GRAVITY),
    }
    for slider in sliders.values():
        slider.on_changed(update_sliders)

    exit_button_ax = plt.axes([0.25, 0, 0.1, 0.03])

    exit_button = Button(exit_button_ax, "Exit")
    exit_button.on_clicked(exit)

    start_time = time.time() - TIME_STEP

    # Simulation loop
    for step in range(NUM_STEPS):
        end_time = time.time()
        TIME_STEP = end_time - start_time
        start_time = end_time
        
        particles, indices, hashes = hash_and_sort(particles)
        densities = compute_density(particles, indices, hashes, SMOOTHING_LENGTH)
        forces = compute_forces(particles, indices, hashes, densities, velocities)

        # Update velocities and positions
        velocities += TIME_STEP * forces / densities[:, None]
        velocities *= 0.8
        particles += TIME_STEP * velocities

        # Enforce boundaries
        enforce_boundaries(particles, velocities, DOMAIN_SIZE)

        # Update visualization
        scatter.set_offsets(particles)
        scatter.set_array(densities)
        ax.set_xlim(0, DOMAIN_SIZE)
        ax.set_ylim(0, DOMAIN_SIZE)
        plt.draw()
        plt.pause(0.0001)

    plt.show()

if __name__ == "__main__":
    with multiprocessing.Pool() as pool:
        main()  # Your main function
