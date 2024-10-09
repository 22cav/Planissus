import random as rdm
import tkinter as tk
from tkinter import Toplevel
import matplotlib.pyplot as plt
import numpy as np
import ttkbootstrap as ttk
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure


# GUI related functions

def update_cursor_label(event):
    """
    This function is linked with the cursor's movement detected by the animated grid, every time that a movement is
    detected, if the animation is still running it just updates the X and Y labels at the bottom of the window. If the
    animation is on pause, then this function also highlight the cell in which the cursor is placed.
    It also manages the case in which the cursor is detected but it is outside of the Grid, printing a message.
    :param event: (MouseEvent) The mouse event triggering the function.
    :return: None
    """

    if not anim_running:
        x, y = event.xdata, event.ydata
        if x is not None and y is not None:
            cursor_label.config(text=f"X: {int(x)}, Y: {int(y)}")
            # Determine the cell index and change its color
            i, j = int(y), int(x)
            grid = Generate_Grid(world, mappa)
            grid[i, j] = 5
            mat.set_data(grid)
            canvas.draw()
        else:
            cursor_label.config(text="Cursor outside the grid")
    else:
        x, y = event.xdata, event.ydata
        if x is not None and y is not None:
            cursor_label.config(text=f"X: {int(x)}, Y: {int(y)}")


def press_pause_button():
    """
    This function is linked to the pressing of the pause button placed in the window.
    When clicked it updates the GUI and it effectively stops the animation showed
    :return: None
    """
    global anim_running, Pause_Button, Restart_Button, New_Map_Button
    if anim_running:
        main_animation.event_source.stop()
        anim_running = False
        Pause_Button.configure(text="Play")
        Restart_Button.configure(state="normal")
        New_Map_Button.configure(state="normal")
    else:
        main_animation.event_source.start()
        anim_running = True
        Pause_Button.configure(text="Stop")
        Restart_Button.configure(state="disabled")
        New_Map_Button.configure(state="disabled")


def restart_sim_button():
    """
    This function is linked to the pressing of the button that lets you restart the simulation.
    The map is the same as the previous but the initialization of the creatures is different.
    The global variables are of great use in this function.
    It also updates the GUI so the buttons are coordinated
    :return: None
    """
    global anim_running, Pause_Button, world, Days_Passed, mappa
    if not anim_running:
        world = World_Generation()
        mat.set_data(Generate_Grid(world, mappa))
        World_Initialization()
        main_animation.event_source.start()
        anim_running = True
        Pause_Button.configure(text="Stop")
        Restart_Button.configure(state="disabled")
        New_Map_Button.configure(state="disabled")
        Days_Passed = 0
    else:
        Guide_Label.configure(text="Pause the simulation in order to restart the world with the same map!")


def inspectCell(event):
    """
    The function is related to the cursor's click when the cursor is placed in the map.
    The function only activates if the cursor is on the map when the animation is on pause and when the cursor is on a
    ground cell, so it doesn't activate if the cursor is on a water cell.
    It detects the x_coordinate and y_coordinate from the map and it calls the graph selection class.
    :param event: The mouse event triggering the function.
    :return: None
    """
    global anim_running, world, mappa
    if event.button == 1 and not anim_running:
        x, y = int(event.xdata), int(event.ydata)
        if x and y:
            if Ground_Check((x, y), mappa):
                graph_selection(x, y, world)


def new_map_button():
    """
    This function is very similar to the restart_sim_button(), but it also creates a new map by calling the same
    function that is used to create a map at the beginning
    :return: None
    """
    global anim_running, Pause_Button, world, mappa, Days_Passed
    if not anim_running:
        mappa = Map_generator()
        world = World_Generation()
        mat.set_data(Generate_Grid(world, mappa))
        World_Initialization()
        main_animation.event_source.start()
        anim_running = True
        Pause_Button.configure(text="Stop")
        Restart_Button.configure(state="disabled")
        New_Map_Button.configure(state="disabled")
        Days_Passed = 0
    else:
        Guide_Label.configure(text="Pause the simulation in order to get a fresh new map!")


def get_vegetob_data(vegetob_list):
    """
    This function is useful only when the user wants to explicitly visualize the data contained in a cell.
    In order to reach this function you must pass through the graph_selection window and select the Visualize_data
    button. It just goes through the list and returns a list with all the interesting values.
    In order to use it properly you have to use as a parameter the list of the creatures you're interested in
    :param vegetob_list: List of Vegetob creatures.
    :return: A string containing data for Vegetob creatures.
    """
    data = ''
    for veg in vegetob_list:
        if not isinstance(veg, Vegetob):
            return ''
        data += '    Vegetob:\n'
        data += f'        density: {veg.get_density()}\n'
        data += f'        days passed since growing: {veg.get_days_passed_since_growing()}\n'
    return data


def get_erbast_data(herd_list):
    """
    same as the get_vegetob_data function, but here the data concerns Erbast creatures
    :param herd_list: List of Herd groups containing Erbast creatures.
    :return: A string containing data for Erbast creatures.
    """
    data = ''
    c = 1
    for herd in herd_list:
        if not isinstance(herd, Herd):
            return ''
        data += f'    Herd {c}:\n'
        k = 1
        for erb in herd.get_members_list():
            if not isinstance(erb, Erbast):
                return ''
            data += f'        Erbast {k}:\n'
            data += f'            energy: {erb.get_energy()}\n'
            data += f'            lifetime: {erb.lifetime}\n'
            data += f'            age: {erb.age}\n'
            data += f'            social attitude: {erb.get_sa()}\n'
            k += 1
        c += 1
    return data


def get_carvix_data(carvix_list):
    """
    same as the get_vegetob_data function, but here the data concerns Carvix creatures
    :param carvix_list: List of Pride groups containing Carvix creatures.
    :return: A string containing data for Carvix creatures.
    """
    data = ''
    c = 1
    for pride in carvix_list:
        if not isinstance(pride, Pride):
            return ''
        data += f'    Pride {c}:\n'
        k = 1
        for car in pride.get_members_list():
            if not isinstance(car, Carvix):
                return ''
            data += f'        Carvix {k}:\n'
            data += f'            energy: {car.get_energy()}\n'
            data += f'            lifetime: {car.lifetime}\n'
            data += f'            age: {car.age}\n'
            data += f'            social attitude: {car.get_sa()}\n'
            k += 1
        c += 1
    return data


def get_cell_data(cell):
    """
    When this function is given a cell contained in world, then it returns a string with all the interesting data
    contained in it. It just performs a string concatenation.
    :param cell: A dictionary representing a cell in the world.
    :return: A string containing data for creatures and elements in the cell.
    """
    data = ''
    if cell.get(Vegetob.__name__):
        data += 'Vegetob list:\n'
        data += get_vegetob_data(cell.get(Vegetob.__name__))
        data += '\n'
    if cell.get(Erbast.__name__):
        data += 'Erbast list:\n'
        data += get_erbast_data(cell.get(Erbast.__name__))
        data += '\n'
    if cell.get(Carvix.__name__):
        data += 'Carvix list:\n'
        data += get_carvix_data(cell.get(Carvix.__name__))
        data += '\n'
    return data


def update(frame):
    """
    This is one of the main function of the whole project.
    It is used to update the function that is shown in the main window, that represents the map, the state of the
    world and of its creatures at a certain instant of time.
    The world lives one day every time this function is called, thanks to the Live_One_Day function.
    The global variables world and mappa are used in order to generate a grid, thanks to the Generate_Grid function.
    The grid is then set as data in the matrix that is shown in the main window.
    :param frame: The frame number for animation.
    :return: A single-member list containing the updated grid.
    """
    global world, mappa
    world = Live_One_Day(world, mappa)
    grid = Generate_Grid(world, mappa)
    mat.set_data(grid)
    return [mat]


# Planissus related functions

def calculate_distance(pos1, pos2):
    """
    It just calculates and return the Manhattan distance between a cell and another based only on their position.
    pos1 and pos2 are the positions of the cells.
    :param pos1: The position of the first cell as a tuple (x, y).
    :param pos2: The position of the second cell as a tuple (x, y).
    :return: The Manhattan distance between the two cells.
    """
    if len(pos1) != len(pos2):
        return None  # Return None for invalid input
    return sum(abs(x - y) for x, y in zip(pos1, pos2))


def calculate_energy_expenditure(distance, group):
    """
    this function is used to calculate how much energy would be used by the whole creature group
    trying to reach a certain cell. The distance is still the Manhattan distance.
    :param distance: The Manhattan distance the group is moving.
    :param group: The group of creatures.
    :return: The calculated energy expenditure.
    """
    return 1 + (distance * group.get_number_of_members())


def calculate_energy_available(group):
    """
    The function is calculating the energy available by the whole group
    :param group: The group of creatures.
    :return: The total available energy for the group.
    """
    return 1 + group.get_total_energy()


def calculate_dangerousness(position):
    """
    The function calculate the dangerousness of a cell based on the presence of Carvix creatures. It is used both in
    Erbast and Carvix because the Carvix creatures could be a danger for both species.
    (Erbasts are hunted and Carvixes can fight each other)
    :param position: The position of the cell as a tuple (x, y).
    :return: The calculated dangerousness level of the cell.
    """
    global world
    danger_level = 1
    if world[position][Carvix.__name__]:
        for pride in world[position][Carvix.__name__]:
            if pride:
                for predator in pride.get_members_list():
                    if not isinstance(predator, Carvix):
                        break
                    danger_level += 1
    return 1 + danger_level


def calculate_conditions_erbast(position):
    """
    This function is used to calculate the amount of preys in a cell, based only on its position.
    It is refered to Erbasts, so the preys in this case are the Vegetobs.
    :param position: The position of the cell as a tuple (x, y).
    :return: The calculated environmental conditions for Erbast in the cell.
    """
    global world
    veg = 0
    if world[position][Vegetob.__name__]:
        for member in world[position][Vegetob.__name__]:
            if not isinstance(member, Vegetob):
                break
            veg += member.get_density()
    return 1 + veg


def calculate_conditions_carvix(position):
    """
    This function is used to calculate the amount of preys in a cell, based only on its position.
    It is refered to Carvixes, so the preys in this case are the Erbasts.
    :param position: The position of the cell as a tuple (x, y).
    :return: The calculated environmental conditions for Carvix in the cell.
    """
    global world
    preys = 0
    if world[position][Erbast.__name__]:
        for member in world[position][Erbast.__name__]:
            if not isinstance(member, Herd):
                break
            preys += member.get_number_of_members()
    return 1 + preys


def Initialize_Vegetob(position):
    """
    One of the first functions of the program.
    It initializes a Vegetob in a given position, randomizing all the other class parameters with some limits.
    It is mainly used in the World Initialization process.
    :param position: The position where the Vegetob is initialized as a tuple (x, y).
    :return: An initialized Vegetob creature.
    """
    density = np.random.randint(0, 20)
    creature = Vegetob(position, density)
    return creature


def Initialize_Erbast(position, energy):
    """
    One of the first functions of the program.
    It initializes a Erbast in a given position, randomizing all the other class parameters with some limits.
    It is mainly used in the World Initialization process.
    :param position: The position where the Erbast is initialized as a tuple (x, y).
    :param energy: The initial energy of the Erbast creature.
    :return: An initialized Erbast creature.
    """
    if not energy:
        energy = rdm.randint(25, int(MAX_ENERGY_ERBAST))
    lifetime = np.random.randint(25, MAX_LIFETIME_ERBAST)
    age = 0
    sa = rdm.randint(1, 5)
    creature = Erbast(position, energy, lifetime, age, sa)
    return creature


def Initialize_Carvix(position, energy):
    """
    One of the first functions of the program.
    It initializes a Carvix in a given position, randomizing all the other class parameters with some limits.
    It is mainly used in the World Initialization process.
    :param position: The position where the Carvix is initialized as a tuple (x, y).
    :param energy: The initial energy of the Carvix creature.
    :return: An initialized Carvix creature.
    """
    if not energy:
        energy = rdm.randint(50, int(MAX_ENERGY_CARVIX))
    lifetime = np.random.randint(25, MAX_LIFETIME_CARVIX)
    age = 0
    sa = rdm.randint(1, 5)
    creature = Carvix(position, energy, lifetime, age, sa)
    return creature


def Random_Position(x_size, y_size):
    """
    As the name suggests, it returns a random position composed of the x and y coordinates.
    The position is composed like: { (x,y) | 1 <= x <= x_size and 1 <= y <= y_size}
    :param x_size: The size of the world in the x-coordinate.
    :param y_size: The size of the world in the y-coordinate.
    :return: A random position as a tuple (x, y) within the world boundaries.
    """
    p = (np.random.randint(1, x_size), np.random.randint(1, y_size))
    return p


def Initialize_Creature(position, creature_type):
    """
    This function is used only in the World Initialization process. It initializes and return a creature having only its
    position and the kind of creature it is.
    :param position: The position where the creature is initialized as a tuple (x, y).
    :param creature_type: The type of creature to initialize (Vegetob, Erbast, or Carvix).
    :return: An initialized creature of the specified type.
    """
    if creature_type == Vegetob:
        creature = Initialize_Vegetob(position)

    elif creature_type == Erbast:
        creature = Initialize_Erbast(position, None)

    elif creature_type == Carvix:
        creature = Initialize_Carvix(position, None)
    else:
        creature = None
    return creature


# In this Planissus world the creatures that move are only organized in groups. The groups can be composed of only
# one member, but it is still classified as one group.

def Initialize_Herd(position, creature):
    """
    In this function a given Erbast creature is initialized in a Herd composed only by a member
    :param position: The position where the Herd is initialized as a tuple (x, y).
    :param creature: The Erbast creature to be added to the Herd.
    :return: An initialized Herd containing the specified Erbast creature.

    """
    herd = Herd(position)
    herd.add_erbast(creature)
    return herd


def Initialize_Pride(position, creature):
    """
    In this function a given Carvix creature is initialized in a Pride composed only by a member
    :param position: The position where the Pride is initialized as a tuple (x, y).
    :param creature: The Carvix creature to be added to the Pride.
    :return: An initialized Pride containing the specified Carvix creature.
    """
    pride = Pride(position)
    pride.add_carvix(creature)
    return pride


def check_list(group_list, class_name):
    """
    This functions works as a check.
    It analyses a group_list or a list of vegetob (Vegetobs are considered as single creatures)
    If the conditions are met, then nothing is changed.
    If the conditions are not met, then the exception gets removed.
    It is a way to ensure that the program works fine.
    :param group_list: The list of creatures in the group.
    :param class_name: The name of the creature class (Vegetob, Erbast, or Carvix).
    :return: None
    """
    if group_list:
        if class_name == Vegetob.__name__:
            for element in group_list:
                if isinstance(element, Vegetob):
                    if element.get_density() <= 0:
                        group_list.remove(element)
                    if element.get_density() > MAX_DENSITY_VEGETOB:
                        element.density = 100
            if group_list:
                while not len(group_list) == 1:
                    if not group_list[1]:
                        break
                    group_list[0].add_density(group_list[1].get_density())
                    group_list.remove(group_list[1])
        elif class_name == Erbast.__name__:
            for element in group_list:
                if isinstance(element, Herd):
                    element.check_members(group_list)
                    if element.check_if_empty():
                        group_list.remove(element)
                    else:
                        for c in element.get_members_list():
                            c.check_conditions()
                else:
                    group_list.remove(element)
        elif class_name == Carvix.__name__:
            for element in group_list:
                if isinstance(element, Pride):
                    element.check_members(group_list)
                    if element.check_if_empty():
                        group_list.remove(element)
                    else:
                        for c in element.get_members_list():
                            c.check_conditions()
                else:
                    group_list.remove(element)
        return


def check_cell(cell):
    """
    Given a cell, it takes every group_list (or list of vegetobs) and it checks it calling the check_list function.
    :param cell: The cell containing group lists of creatures.
    :return: None
    """
    check_list(cell[Vegetob.__name__], Vegetob.__name__)
    check_list(cell[Erbast.__name__], Erbast.__name__)
    check_list(cell[Carvix.__name__], Carvix.__name__)
    return


def track_function(position):
    """
    This function can be ignored.
    It is not present in the overall program as it has some limitations.
    Its initial purpose was to track a group of moving creatures and show its trajectory on a map where only the group's
    movements were shown
    :param position: The position to track as a tuple (x, y).
    :return: None
    """
    p = rdm.choice(world[position][Erbast.__name__] + world[position][Carvix.__name__])
    p.track = True

    fig2, ax2 = plt.subplots()
    ax2.axis('off')

    grid = np.zeros(world_size)
    grid[mirror(p.position[0], p.position[1])] = 1
    data2 = [p.position]
    mat2 = ax2.matshow(grid, vmin=0, vmax=1)

    def update_traj(frame):
        grid[mirror(p.position[0], p.position[1])] = 1
        if len(data2) > 5:
            data2.remove(data2[0])
        data2.append(p.position)
        mat2.set_data(grid)
        plt.show()
        return [mat2]

    traj_anim = FuncAnimation(fig2, update_traj, frames=100, interval=400)
    plt.show()


def Generate_Grid(world, mappa):
    """
    This function is very important to the program as it generates the grid that is shown effectively on the screen.
    Its working is not very complex, as it assigns a value between 0 and 4 to every position based on the presence of
    creatures in the world. It also detects whether there is ground or water in the map.
    :param world: The world containing creatures and their positions.
    :param mappa: The map representing environmental features.
    :return: A grid (numpy array) representation of the world with values indicating creature presence and features.
    """
    grid = np.zeros(world_size)
    for x in range(world_size[0]):
        for y in range(world_size[1]):
            if not Ground_Check((x, y), mappa):
                grid[mirror(x, y)] = 0
            else:
                if world[x, y][Carvix.__name__]:
                    if isinstance(world[x, y][Carvix.__name__][0], Pride):
                        grid[mirror(x, y)] = 4
                elif world[x, y][Erbast.__name__]:
                    if isinstance(world[x, y][Erbast.__name__][0], Herd):
                        grid[mirror(x, y)] = 3
                elif world[x, y][Vegetob.__name__]:
                    if isinstance(world[x, y][Vegetob.__name__][0], Vegetob):
                        grid[mirror(x, y)] = 2
                else:
                    grid[mirror(x, y)] = 1
    return grid


def World_Generation():
    """
    This function just generates an empty "world" variable. Nothing is put inside of it yet.
    :return: An empty world structure.
    """
    global world_size
    _world = {(x, y): {'Vegetob': [], 'Erbast': [], 'Carvix': []} for x in range(world_size[0]) for y in
              range(world_size[1])}
    return _world


def World_Initialization():
    """
    Initializes the world with creatures.
    This function populates the world with creatures, considering their types and positions on the map.
    :return: None
    """
    global num_creatures, world, mappa

    for _ in range(num_creatures):

        # Select a random creature type
        if rdm.random() < 0.5:
            creature_type = creature_types[0]
        else:
            if rdm.random() < 0.8:
                creature_type = creature_types[1]
            else:
                creature_type = creature_types[2]

        # Select some random Coordinates
        position = Random_Position(world_size[0], world_size[1])

        # Check if the position is suitable (on ground)
        while not Ground_Check(position, mappa):
            position = Random_Position(world_size[0], world_size[1])

        # Initialize a creature
        creature = Initialize_Creature(position, creature_type)
        cell = world[position]
        creature_list = cell[creature_type.__name__]

        # Initialize the group
        if not creature_list:
            if creature_type == Erbast:
                herd = Initialize_Herd(position, creature)
                creature_list.append(herd)
            elif creature_type == Carvix:
                pride = Initialize_Pride(position, creature)
                creature_list.append(pride)
            else:
                creature_list.append(creature)
        else:
            if creature_type == Erbast:
                creature_list[0].add_erbast(creature)
            elif creature_type == Carvix:
                creature_list[0].add_carvix(creature)
            else:
                creature_list[0].add_density(1)
    return


def mirror(x, y):
    return y, x


def Map_generator():
    """
    Generates a map with water and ground cells.
    This function generates a map with water ('W') and ground ('G') cells. It ensures that the borders are water cells
    and then randomly adds more water cells to the map.
    :return: A numpy array representing the generated map.
    """
    global world_size, total_cells

    # Create an empty map filled with ground ('G')
    _mappa = np.full(world_size, 'G', dtype=str)

    # Define the border as a slice and set border cells to 'W'
    border = np.zeros(world_size, dtype=bool)
    border[0, :] = True
    border[-1, :] = True
    border[:, 0] = True
    border[:, -1] = True
    _mappa[border] = 'W'

    # Get the coordinates of water cells in the border
    water_cells = list(zip(*np.where(border)))

    # Determine the number of water cells to add
    num_water_cells = np.random.randint(total_cells * 0.4, total_cells * 0.5)

    def is_within_bounds(x, y):
        return 0 <= x < world_size[0] and 0 <= y < world_size[1]

    # Add water cells until the desired count is reached
    while not len(water_cells) == num_water_cells:
        possible_additions = []
        while not possible_additions:
            w = rdm.choice(water_cells)
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    x, y = w[0] + dx, w[1] + dy
                    if is_within_bounds(x, y) and (x, y) not in water_cells:
                        possible_additions.append((x, y))
        r = rdm.choice(possible_additions)
        water_cells.append(r)
        _mappa[r] = "W"

    # Remove single-cell islands
    _mappa = Delete_Single_Cell_Islands(_mappa)

    return _mappa


def Delete_Single_Cell_Islands(mappa):
    """
    Deletes single-cell islands from the map.
    This function checks each cell in the map and if it's a ground cell ('G') and it's completely surrounded by
    other ground cells, it converts it to a water cell ('W').
    :param mappa: The map to process.
    :return: The modified map with single-cell islands removed.
    """
    for i in range(world_size[0]):
        for j in range(world_size[1]):
            if mappa[i, j] == 'G' and 'G' not in [mappa[i + 1, j], mappa[i - 1, j], mappa[i, j + 1], mappa[i, j - 1]]:
                # If a ground cell doesn't have at least one ground neighbour, change it to water
                mappa[i, j] = 'W'
    return mappa


def Random_Growth(world, mappa):
    """
    Simulates random growth of vegetation on the world map. This function attempts to grow vegetation at random
    positions on the map. It checks if the position is suitable for growth and increases the density of existing
    vegetation or initializes new vegetation if none exists at that position. :param world: The world data structure
    containing information about creatures and vegetation.
    :param mappa: The map that defines suitable locations for growth.
    :return: The modified world data structure after simulating random growth.
    """
    global TotGrowth

    # Perform random growth for a specified number of repetitions
    for rep in range(TotGrowth):
        while True:
            # Select a random position for growth
            p = Random_Position(world_size[0], world_size[1])
            if Ground_Check(p, mappa):
                break

        # Check if there is existing vegetation at the selected position
        if len(world[p][Vegetob.__name__]) == 0:
            # If no vegetation exists, initialize new vegetation
            world[p][Vegetob.__name__] = [Initialize_Vegetob(p)]
        else:
            # If vegetation already exists, increase its density
            creature = world[p][Vegetob.__name__][0]
            creature.density += 2

    return world


def Ground_Check(position, mappa):
    """
    Checks if a given position on the map is on the ground.
    This function takes a position and checks whether the corresponding cell on the map ('mappa')
    is 'G' (ground) or not.
    :param position: The position to check on the map.
    :param mappa: The map that defines ground and non-ground cells.
    :return: True if the position is on the ground, False otherwise.
    """
    if mappa[position] == 'G':
        return True
    else:
        return False


def herd_union(herd2, group_list):
    """
    Performs a union operation between two herds in a group list.
    This function attempts to merge two herds ('herd1' and 'herd2') in the 'group_list' if certain conditions are met.
    If the conditions are not met, it performs various operations such as removing members or creating new herds.
    :param herd2: The second herd to be considered for union.
    :param group_list: The list containing herds or groups of creatures.
    :return: None
    """
    # Get the first herd in the group list
    herd1 = group_list[0]

    # Check if 'herd1' and 'herd2' are instances of the 'Herd' class and share the same position
    if not (isinstance(herd1, Herd) and isinstance(herd2, Herd)) \
            or not herd1.position == herd2.position:
        # If conditions are not met, perform checks and return
        check_list(herd1, Erbast.__name__)
        check_list(herd2, Erbast.__name__)
        return

    # Attempt to merge the members of 'herd2' into 'herd1' if the member count in 'herd1' is less than 15
    for e in herd2.get_members_list():
        if len(herd1.get_members_list()) < 15:
            herd2.remove_member(e)
            herd1.add_erbast(e)
        else:
            # If 'herd1' is already full, remove and create a new herd for the member
            herd2.remove_and_create_herd(e, group_list)

    # Remove 'herd2' from the group list
    group_list.remove(herd2)

    # Perform checks on the group list
    check_list(group_list, Erbast.__name__)

    return


def pride_union(pride2, group_list):
    """
    Performs a union operation between two prides in a group list.
    This function attempts to merge two prides ('pride1' and 'pride2') in the 'group_list' if certain conditions are
    met. If the conditions are not met, it performs various operations such as removing members or creating new prides.
    :param pride2: The second pride to be considered for union.
    :param group_list: The list containing prides or groups of creatures.
    :return: None
    """
    # Get the first pride in the group list
    pride1 = group_list[0]

    # Check if 'pride1' and 'pride2' are instances of the 'Pride' class and share the same position
    if not (isinstance(pride1, Pride) and isinstance(pride2, Pride)) \
            or not pride1.position == pride2.position:
        # If conditions are not met, perform checks and return
        check_list(pride1, Carvix.__name__)
        check_list(pride2, Carvix.__name__)
        return

    # Attempt to merge the members of 'pride2' into 'pride1' if the member count in 'pride1' is less than 7
    for e in pride2.get_members_list():
        if len(pride1.get_members_list()) < 7:
            pride2.remove_member(e)
            pride1.add_carvix(e)
        else:
            # If 'pride1' is already full, remove and create a new pride for the member
            pride2.remove_and_create_pride(e, group_list)

    # Remove 'pride2' from the group list
    group_list.remove(pride2)

    # Perform checks on the group list
    check_list(group_list, Carvix.__name__)

    return


def opposite_direction(direction):
    """
    Calculates the opposite direction for a given direction.
    This function takes a direction as a tuple (dx, dy) and returns the opposite direction.
    :param direction: The direction for which to find the opposite.
    :return: The opposite direction as a tuple (dx, dy).
    """
    if direction == (0, 1):
        return 0, -1
    elif direction == (0, -1):
        return 0, 1
    elif direction == (1, 0):
        return -1, 0
    elif direction == (-1, 0):
        return 1, 0
    else:
        return None


def pride_fight(pride1, pride2, pride_list):
    """
    Simulates a fight between two prides of creatures.
    This function simulates a fight between two prides ('pride1' and 'pride2') in a list of prides ('pride_list').
    It checks the energy levels of the leaders of both prides and performs actions accordingly, including removing
    members and prides from the list.
    :param pride1: The first pride involved in the fight.
    :param pride2: The second pride involved in the fight.
    :param pride_list: The list containing all prides of creatures.
    :return: The surviving pride or None if both prides are eliminated.
    """
    # Check if 'pride1' and 'pride2' are instances of the 'Pride' class
    if not (isinstance(pride1, Pride) and isinstance(pride2, Pride)):
        return

    # Check and update members of 'pride1' and 'pride2'
    pride1.check_members(pride_list)
    pride2.check_members(pride_list)

    # Continue the fight as long as both prides have leaders and members
    while pride1 and pride2 and pride1.get_members_list() and pride2.get_members_list():
        l1 = pride1.get_leader()
        l2 = pride2.get_leader()

        # Compare the energy levels of the leaders
        if l1.get_energy() > l2.get_energy():
            l1.rem_energy(l2.get_energy())
            for member in pride2.get_members_list():
                pride2.remove_member(member)
            return pride1
        elif l1.get_energy() < l2.get_energy():
            l2.rem_energy(l1.get_energy())
            for member in pride1.get_members_list():
                pride1.remove_member(member)
            return pride2
        else:
            # If leaders have equal energy, remove them from their respective prides
            pride1.remove_member(l1)
            pride2.remove_member(l2)
            pride1.check_members(pride_list)
            pride2.check_members(pride_list)

    # Check the outcome of the fight and remove eliminated prides
    if not pride1:
        world[pride1.position][Carvix.__name__].remove(pride1)
        return pride2
    elif not pride2:
        world[pride2.position][Carvix.__name__].remove(pride2)
        return pride1
    else:
        return None


def prides_HOS(position, pride_list, world, mappa):
    """
    Handles the hunting and merging of prides at a given position.
    This function manages the behavior of prides at a specific position, including hunting for prey, merging prides if
    conditions are met, and handling fights between prides.
    :param position: The position on the map where prides are interacting.
    :param pride_list: The list of prides in the vicinity.
    :param world: The world data structure containing information about creatures and vegetation.
    :param mappa: The map that defines suitable locations for movement and interactions.
    :return: None
    """
    # Check the integrity of the pride list
    check_list(pride_list, Carvix.__name__)

    # If there are no prides, return
    if not pride_list:
        return

    # If there is only one pride in the list, handle its behavior
    if len(pride_list) == 1:
        hunters = world[position][Carvix.__name__][0]
        targets = world[position][Erbast.__name__]

        # Check if the sole pride in the list is an instance of 'Pride'
        if isinstance(hunters, Pride):
            if not targets:
                # If no targets are present, move the pride
                hunters.move_pride(world, mappa)
            else:
                # If targets are present, initiate hunting
                hunters.hunt(targets)

    else:
        # Continue handling prides as long as there is more than one pride in the list
        while not len(pride_list) == 1:
            if not pride_list:
                return

            # Randomly select a pride from the list
            p = rdm.choice(pride_list)

            # Determine whether to merge prides or initiate a fight based on a probability check
            if rdm.random() > p.get_total_energy() / 300:
                # Merge prides if the probability condition is met
                pride_to_be_removed = rdm.choice(pride_list[1:])
                pride_union(pride_to_be_removed, pride_list)
            else:
                # Choose two prides for a fight and remove them from the list
                p1 = p
                pride_list.remove(p1)
                p2 = rdm.choice(pride_list)
                pride_list.remove(p2)

                # Initiate a pride fight and determine the winner
                winner = pride_fight(p1, p2, world[position][Carvix.__name__])

                # Add the winner (if any) back to the pride list
                if winner:
                    pride_list.append(winner)

        hunters = world[position][Carvix.__name__][0]
        targets = world[position][Erbast.__name__]

        # Check if the sole remaining pride in the list is an instance of 'Pride'
        if isinstance(hunters, Pride):
            if not targets:
                # If no targets are present, move the pride
                hunters.move_pride(world, mappa)
            else:
                # If targets are present, initiate hunting
                hunters.hunt(targets)

    return


def Live_One_Day(world, mappa):
    """
    Simulates one day of activity in the world.
    This function simulates various activities for creatures and vegetation in the world, including growth, movement,
    grazing, struggling, spawning, and cell checks.
    :param world: The world data structure containing information about creatures and vegetation.
    :param mappa: The map that defines suitable locations for movement and interactions.
    :return: The updated world data structure after simulating one day.
    """
    # GROWING
    # Iterate through the world
    for i in range(world_size[0]):
        for j in range(world_size[1]):

            # If there is no ground continue the iteration
            if not Ground_Check((i, j), mappa):
                continue

            # If there exist a vegetob in the cell, then call its method grow
            if world[i, j][Vegetob.__name__]:
                if isinstance(world[i, j][Vegetob.__name__][0], Vegetob):
                    world[i, j][Vegetob.__name__][0].grow()

    # Call the Random_Growth function that let Vegetobs grow from nothing in the world
    world = Random_Growth(world, mappa)

    # MOVEMENT
    # Iterate through the world
    for i in range(world_size[0]):
        for j in range(world_size[1]):

            # If there is no ground continue the iteration
            if not Ground_Check((i, j), mappa):
                continue

            # If there is at least one Herd in the cell then decide whether to move it or not
            if world[i, j][Erbast.__name__]:
                while True:

                    # All the white herds contained in the cell are stored in a list
                    white_herds = [herd for herd in world[i, j][Erbast.__name__] if herd.get_color() == "White"]

                    # Interrupt the loop when the white herds are finished
                    if not white_herds:
                        break

                    # The white herds can either move (they become Black) or stay still (they become Green)
                    for herd in white_herds:
                        if isinstance(herd, Herd):
                            herd.social_decision(world, mappa)
                            herd.check_members(world[i, j][Erbast.__name__])

                        # If the herd is empty, then it is set black and kept waiting until the next check
                        if herd.check_if_empty():
                            herd.set_color("Black")

            # If there is at least one Pride in the cell then decide whether to move it or not
            if world[i, j][Carvix.__name__]:
                while True:

                    # All the white prides contained in the cell are stored in a list
                    white_prides = [pride for pride in world[i, j][Carvix.__name__] if pride.get_color() == "White"]

                    # Interrupt the loop when the white prides are finished
                    if not white_prides:
                        break

                    # The white prides can either move (they become Black) or stay still (they become Green)
                    for pride in white_prides:
                        if isinstance(pride, Pride):
                            pride.social_decision(world, mappa)
                            pride.check_members(world[i, j][Carvix.__name__])

                        # If the pride is empty, then it is set black and kept waiting until the next check
                        if pride.check_if_empty():
                            pride.set_color("Black")

    # GRAZING AND STRUGGLING
    # Iterate through the world
    for i in range(world_size[0]):
        for j in range(world_size[1]):

            # Continue the iteration if it is a Water cell
            if not Ground_Check((i, j), mappa):
                continue

            # Check the cell
            check_cell(world[i, j])

            # If there is at least one Herd
            if world[i, j][Erbast.__name__]:

                # The green herds are stored in a list (The ones that did not move before)
                green_herds = [herd for herd in world[i, j][Erbast.__name__] if herd.get_color() == "Green"]

                # The green Herds now can graze
                if green_herds:
                    for herd in green_herds:
                        check_list(green_herds, Erbast.__name__)
                        herd.graze(world)
                        herd.set_color("Black")
                        herd.check_members(world[i, j][Erbast.__name__])

                    # If there is more than one herd, then unite every herd in one
                    if len(world[i, j][Erbast.__name__]) > 1:
                        for herd_to_be_removed in world[i, j][Erbast.__name__][1:]:
                            herd_union(herd_to_be_removed, world[i, j][Erbast.__name__])

    # Iterate through the world
    for i in range(world_size[0]):
        for j in range(world_size[1]):

            # Continue the iteration if it is a Water cell
            if not Ground_Check((i, j), mappa):
                continue

            # If there is at least one Pride
            if world[i, j][Carvix.__name__]:

                # Store the green prides in a list (The prides that did not move before)
                green_prides = [pride for pride in world[i, j][Carvix.__name__] if pride.get_color() == "Green"]
                if green_prides:
                    check_list(green_prides, Carvix.__name__)

                    # The prides can Fight, Hunt, Unite
                    prides_HOS((i, j), green_prides, world, mappa)

    # SPAWNING
    # Iterate through the world
    for i in range(world_size[0]):
        for j in range(world_size[1]):

            # Continue the iteration if it is a Water cell
            if not Ground_Check((i, j), mappa):
                continue

            # If there is at least one Herd
            if world[i, j][Erbast.__name__]:

                # Iterates for every Herd
                for erbast_instance in world[i, j][Erbast.__name__]:
                    if isinstance(erbast_instance, Herd):
                        # Let the Herd spawn the children and add age to its members
                        erbast_instance.Spawn()

                        # Check the members and their stats
                        erbast_instance.check_members(world[i, j][Erbast.__name__])

                        # Set the color of the herd back to white
                        erbast_instance.reset_color()

            # If there is at least one Pride
            if world[i, j][Carvix.__name__]:
                for carvix_instance in world[i, j][Carvix.__name__]:
                    if isinstance(carvix_instance, Pride):
                        # Let the Pride spawn the children and add age to its members
                        carvix_instance.Spawn()

                        # Check the members and their stats
                        carvix_instance.check_members(world[i, j][Carvix.__name__])

                        # Set the color of the pride back to white
                        carvix_instance.reset_color()

    # CHECK EVERY CELL
    # Iterate through the world
    for i in range(world_size[0]):
        for j in range(world_size[1]):
            if Ground_Check((i, j), mappa):
                # Check that every cell is normal
                check_cell(world[i, j])

    global Days_Passed

    # Update the days passed and the label in the main window
    Days_Passed += 1
    days_passed_label.configure(text=f"Days Passed: {Days_Passed}")

    return world


def get_data(_world, _position):
    """
    Get data from a specific cell in the world.
    This function retrieves data from a specific cell in the world and returns it as a list.
    :param _world: The world data structure containing information about creatures and vegetation.
    :param _position: The position (coordinates) of the cell from which to retrieve data.
    :return: A list containing three values:
             - d1: Total density of Vegetob instances in the cell.
             - d2: Total number of creatures in Erbast herds in the cell.
             - d3: Total number of creatures in Carvix prides in the cell.
    """
    # Initialize variables to store data
    d1, d2, d3 = 0, 0, 0

    # Calculate the total density of Vegetob instances in the cell
    for vegetob in _world[_position][Vegetob.__name__]:
        if isinstance(vegetob, Vegetob):
            d1 += vegetob.get_density()

    # Calculate the total number of creatures in Erbast herds in the cell
    for herd in _world[_position][Erbast.__name__]:
        if herd and isinstance(herd, Herd):
            d2 += len(herd.get_members_list())

    # Calculate the total number of creatures in Carvix prides in the cell
    for pride in _world[_position][Carvix.__name__]:
        if pride and isinstance(pride, Pride):
            d3 += len(pride.get_members_list())

    # Return a list containing the calculated data
    return [d1, d2, d3]


# GUI classes

class visualize_cell_data:
    def __init__(self, position, world1):
        # Initialize the visualization object with the current time, cell position, and world data
        self.time = Days_Passed
        self.position = position
        self.cell = world1[self.position]

        # Create the main application window
        window = tk.Tk()
        window.title('Data visualization')

        # Create a title label to display information about the cell and time
        title = tk.Label(window, text=f'This is the data contained in the cell {self.position},'
                                      f' at the instant of time {self.time}')
        title.pack(padx=20, pady=20)

        # Create a text widget to display cell data
        text = tk.Text(window)
        text.pack(pady=20)

        # Retrieve data for the cell using the 'get_cell_data' function
        data = get_cell_data(self.cell)

        # Insert the retrieved data into the text widget
        text.insert(tk.END, data)

        # Start the main GUI event loop
        window.mainloop()


class graph_selection:
    def __init__(self, x, y, world1):
        # Initialize instance variables
        self.labels = ['Vegetob', 'Erbast', 'Carvix']
        self.ani1 = None
        self.data = None
        self.fig1 = None
        self.ax1 = None
        self.world = world1
        self.position = (int(x), int(y))
        self.selection_window = Toplevel()
        self.selection_window.resizable(False, False)
        self.selection_window.geometry('600x400')
        self.selection_window.title('Select graph')

        # Create labels and buttons for graph selection
        self.label = tk.Label(self.selection_window, text='Which type of graph\nwould you like to see?',
                              font=("Fixedsys", 15), fg="white", )
        self.label.pack(padx=20, pady=10)

        # Create a frame to hold the buttons
        self.button_space = tk.Frame(master=self.selection_window)

        """self.track_button = tk.Button(self.button_space, text='Track a Random\ngroup in the Cell',
                                      command = self.press_track, font=('Fixedsys', 10))
        self.track_button.pack(side='top')"""

        # Create a frame to hold the top buttons
        self.top_button_frame = tk.Frame(self.button_space)

        # Create the "Visualize Data" button
        self.visualize_data_button = tk.Button(self.top_button_frame, text='Visualize\nthe Data',
                                               command=self.press_visualize, font=('Fixedsys', 10))
        self.visualize_data_button.grid(row=0, column=0, padx=20, pady=5, sticky='nsew')

        # Create the "Histogram" button
        self.hist_button = tk.Button(self.top_button_frame, text='Real Time\nHistogram', command=self.press_hist,
                                     font=('Fixedsys', 10))
        self.hist_button.grid(row=0, column=1, padx=20, pady=5, sticky='nsew')

        # Create the "Graph" button
        self.graph_button = tk.Button(self.top_button_frame, text='Real Time\nGraph', command=self.press_graph,
                                      font=('Fixedsys', 10))
        self.graph_button.grid(row=0, column=2, padx=20, pady=5, sticky='nsew')

        self.top_button_frame.grid_columnconfigure(0, weight=1)
        self.top_button_frame.grid_columnconfigure(1, weight=1)
        self.top_button_frame.grid_columnconfigure(2, weight=1)
        self.top_button_frame.pack(padx=20, pady=20, fill='x')

        # Create buttons for creature type selection
        self.button2_space = tk.Frame(self.button_space)
        self.vegetob_button = tk.Button(self.button2_space, text='  Vegetob  ', state="disabled",
                                        command=self.press_vegetob, font=('Fixedsys', 10))
        self.vegetob_button.pack(padx=30, pady=10, side='left')
        self.erbast_button = tk.Button(self.button2_space, text='  Erbast  ', state="disabled",
                                       command=self.press_erbast, font=('Fixedsys', 10))
        self.erbast_button.pack(padx=30, pady=10, side='left')
        self.carvix_button = tk.Button(self.button2_space, text='  Carvix  ', state="disabled",
                                       command=self.press_carvix, font=('Fixedsys', 10), padx=10)
        self.carvix_button.pack(padx=30, pady=10, side='left')
        self.button2_space.pack(padx=20, pady=20)

        # Create a confirm button
        self.confirm_button = tk.Button(self.button_space, text='Confirm', state="disabled", command=self.press_confirm,
                                        font=('Fixedsys', 10))
        self.confirm_button.pack(side='bottom', pady=20)

        # Pack the button space
        self.button_space.pack(expand=True)

        # Display the cell position
        self.position_label = tk.Label(self.selection_window, text=f'The position of the cell is: {int(x)}, {int(y)}',
                                       font=("Courier", 10), fg="white", )
        self.position_label.pack(padx=20, side="bottom")
        self.selection_window.mainloop()

    # Define button press functions

    def press_track(self):
        # Ignore, supposed to be a function to track a random group in the cell
        self.selection_window.destroy()
        track_function(self.position)
        return

    def press_visualize(self):
        # Function to visualize cell data
        self.selection_window.destroy()
        visualize_cell_data(self.position, self.world)

    def press_hist(self):
        # Function to prepare for histogram display
        self.hist_button.config(relief=tk.SUNKEN)
        self.graph_button.config(relief=tk.RAISED)
        self.visualize_data_button.config(relief=tk.SUNKEN)
        self.carvix_button.configure(state="disabled")
        self.erbast_button.configure(state="disabled")
        self.vegetob_button.configure(state="disabled")
        self.confirm_button.configure(state="active")

    def press_graph(self):
        # Function to prepare for real-time graph display
        self.hist_button.config(relief=tk.RAISED)
        self.graph_button.config(relief=tk.SUNKEN)
        self.visualize_data_button.config(relief=tk.SUNKEN)
        self.carvix_button.configure(state="normal")
        self.erbast_button.configure(state="normal")
        self.vegetob_button.configure(state="normal")
        self.confirm_button.configure(state="disabled")

    def press_carvix(self):
        # Function to select Carvix for real-time graph
        self.carvix_button.config(relief=tk.SUNKEN)
        self.erbast_button.config(relief=tk.RAISED)
        self.vegetob_button.config(relief=tk.RAISED)
        self.confirm_button.configure(state="active")

    def press_vegetob(self):
        # Function to select Vegetob for real-time graph
        self.carvix_button.config(relief=tk.RAISED)
        self.erbast_button.config(relief=tk.RAISED)
        self.vegetob_button.config(relief=tk.SUNKEN)
        self.confirm_button.configure(state="active")

    def press_erbast(self):
        # Function to select Erbast for real-time graph
        self.carvix_button.config(relief=tk.RAISED)
        self.erbast_button.config(relief=tk.SUNKEN)
        self.vegetob_button.config(relief=tk.RAISED)
        self.confirm_button.configure(state="active")

    def press_confirm(self):
        # Function to confirm and display the selected graph
        self.fig1, self.ax1 = plt.subplots()
        if self.hist_button['relief'] == 'sunken':
            # Display histogram
            self.data = get_data(self.world, self.position)
            self.ax1.bar(self.labels, self.data, color=['#006400', '#FFB90F', '#cc2702'])
            for i, val in enumerate(self.data):
                plt.text(i, 1, str(val), ha='center')
            plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
            plt.xlabel('Data')
            plt.ylabel('Value')
            self.ani1 = FuncAnimation(self.fig1, self.update_hist_graph, interval=400, frames=100)
            self.selection_window.destroy()
            plt.show()
        elif self.graph_button['relief'] == 'sunken':
            # Display real-time graph based on creature type selection
            if self.carvix_button['relief'] == 'sunken':
                self.data = [get_data(self.world, self.position)[2]]
                self.ani1 = FuncAnimation(self.fig1, self.update_c_graph, interval=400, frames=100)
                plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
                plt.xlabel('Time')
                plt.ylabel('Value')
                self.selection_window.destroy()
                plt.show()
            elif self.erbast_button['relief'] == 'sunken':
                self.data = [get_data(self.world, self.position)[1]]
                self.ani1 = FuncAnimation(self.fig1, self.update_e_graph, interval=400, frames=100)
                plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
                plt.xlabel('Time')
                plt.ylabel('Value')
                self.selection_window.destroy()
                plt.show()
            elif self.vegetob_button['relief'] == 'sunken':
                self.data = [get_data(self.world, self.position)[0]]
                self.ani1 = FuncAnimation(self.fig1, self.update_v_graph, interval=400, frames=100)
                plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
                plt.xlabel('Time')
                plt.ylabel('Value')
                self.selection_window.destroy()
                plt.show()
            else:
                pass
        return

    # Define update functions for real-time graphs

    def update_hist_graph(self, frame):
        global world
        self.world = world
        self.data = get_data(self.world, self.position)
        self.ax1.clear()
        self.ax1.bar(self.labels, self.data, color=['#006400', '#FFB90F', '#cc2702'])
        for i, val in enumerate(self.data):
            plt.text(i, 1, str(val), ha='center')
        plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
        plt.xlabel('Data')
        plt.ylabel('Value')
        return self.ax1

    def update_c_graph(self, frame):
        global world, anim_running
        self.world = world
        if anim_running:
            self.data.append(get_data(self.world, self.position)[2])
        self.ax1.clear()
        self.ax1.plot(self.data)
        plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
        plt.xlabel('Time')
        plt.ylabel('Value')
        return self.ax1

    def update_e_graph(self, frame):
        global world, anim_running
        self.world = world
        if anim_running:
            self.data.append(get_data(self.world, self.position)[1])
        self.ax1.clear()
        self.ax1.plot(self.data)
        plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
        plt.xlabel('Time')
        plt.ylabel('Value')
        return self.ax1

    def update_v_graph(self, frame):
        global world, anim_running
        self.world = world
        if anim_running:
            self.data.append(get_data(self.world, self.position)[0])
        self.ax1.clear()
        self.ax1.plot(self.data)
        plt.title(f'Data from the cell ({int(self.position[0])}, {int(self.position[1])})')
        plt.xlabel('Time')
        plt.ylabel('Value')
        return self.ax1


# Planissus related classes


# Define the Herd class
class Herd:
    def __init__(self, position):
        # Initialize a Herd with a target, position, and attributes
        self.target = None  # Initialize target as None
        self.position = position  # Set the initial position
        self.ErbastList = []  # Initialize a list to store Erbast instances
        self.color = "White"  # Initialize the color attribute as "White"
        self.direction = None  # Initialize direction as None
        self.strategy = rdm.choice(Erbast_Strategies + ['None'])  # Randomly select a strategy
        self.track = False  # Initialize track as False

    # Getter method to retrieve the direction
    def get_direction(self):
        return self.direction

    # Setter method to set the direction
    def set_direction(self, d):
        self.direction = d
        return

    # Check if a direction is set
    def check_direction(self):
        if self.direction:
            return True
        else:
            return False

    # Getter method to retrieve the color
    def get_color(self):
        return self.color

    # Setter method to set the color
    def set_color(self, color):
        self.color = color
        return

    # Reset the color to "White"
    def reset_color(self):
        self.color = "White"
        return

    # Add an Erbast instance to the Herd
    def add_erbast(self, erbast):
        if isinstance(erbast, Erbast):
            self.ErbastList.append(erbast)
        return

    # Get the number of Erbast members in the Herd
    def get_number_of_members(self):
        return len(self.ErbastList)

    # Get the list of Erbast members
    def get_members_list(self):
        return self.ErbastList

    # Remove a member from the Herd
    def remove_member(self, member):
        if member in self.ErbastList:
            self.ErbastList.remove(member)

    # Get the leader of the Herd based on energy and sa (social aggression)
    def get_leader(self):
        if self.ErbastList:
            lead = max(self.ErbastList, key=lambda x: x.energy + x.sa)
            self.ErbastList.remove(lead)
            self.ErbastList.insert(0, lead)
            return lead

    # Get the total energy of the Herd by summing the energy of its members
    def get_total_energy(self):
        return sum(member.get_energy() for member in self.get_members_list() or [])

    # Move the Herd in the world based on certain conditions
    def move_herd(self, world, mappa):
        # Store the old position
        old_position = self.position
        directions = [(1, 0), (-1, 0), (0, 1), (0, -1), (0, 0)]  # Possible directions
        weighted_directions = []  # Store valid directions

        # Iterate over possible directions
        for direction in directions:
            new_position = (old_position[0] + direction[0], old_position[1] + direction[1])

            # Check if the new position is within world boundaries
            if not (0 <= new_position[0] < world_size[0] and 0 <= new_position[1] < world_size[1]):
                continue

            # Check if the ground at the new position is suitable
            if Ground_Check(new_position, mappa):
                weighted_directions.append(direction)

        # If there are valid directions, select the one with the highest evaluation
        if weighted_directions:
            selected_direction = max(weighted_directions, key=lambda d: self.evaluate_direction(d))

            # If the selected direction is (0, 0), set the color to 'Green'
            if selected_direction == (0, 0):
                self.set_color('Green')
                return

            # Update the position and move the Herd
            new_position = (old_position[0] + selected_direction[0], old_position[1] + selected_direction[1])
            world[self.position][Erbast.__name__].remove(self)
            self.set_direction(selected_direction)
            self.position = new_position
            world[new_position][Erbast.__name__].append(self)

            # Update the position of the Erbast members in the Herd
            for c in self.ErbastList:
                c.rem_energy(MOVEMENT_COST)
                c.set_position(new_position)

    # Evaluate the desirability of a direction for the Herd's movement
    def evaluate_direction(self, direction):
        weighted_sum = 0
        # Calculate the new position based on the direction
        new_position = (self.position[0] + direction[0], self.position[1] + direction[1])

        # Check if the new position is within world boundaries
        if not (0 <= new_position[0] < world_size[0] and 0 <= new_position[1] < world_size[1]):
            return weighted_sum

        danger_level = 0
        # Iterate over a neighborhood around the new position
        for i in range(self.position[0] + direction[0] - Awareness, self.position[0] + direction[0] + Awareness + 1):
            for j in range(self.position[1] + direction[1] - Awareness,
                           self.position[1] + direction[1] + Awareness + 1):
                if not (0 <= i < world_size[0] and 0 <= j < world_size[1]):
                    continue
                if Ground_Check((i, j), mappa):
                    danger_level += calculate_dangerousness((i, j))
                    if not self.target:
                        continue
                    if calculate_conditions_erbast((i, j)) > calculate_conditions_erbast(self.target):
                        self.target = (i, j)

        # Calculate energy expenditure and target density
        if self.target:
            energy_expenditure = calculate_energy_expenditure(calculate_distance(self.position, self.target), self)
            target_density = calculate_conditions_erbast(self.target)
        else:
            target_density, energy_expenditure = 1, 1

        # Group Behavior: Factor in the proximity of other herd members
        group_energy = calculate_energy_available(self)
        vegetob_density = calculate_conditions_erbast(new_position)

        if direction == self.direction:
            dir_weight = 1
        elif direction == opposite_direction(self.direction):
            dir_weight = 0.25
        else:
            dir_weight = 0.75

        if self.strategy == 'Escape':
            danger_level = danger_level * 3
        elif self.strategy == 'Fight':
            danger_level = danger_level * 0.5

        # Combine all factors with appropriate weights to compute a weighted sum
        weighted_sum = (
                dir_weight * target_density * vegetob_density * group_energy / energy_expenditure * danger_level
        )
        return weighted_sum

    # Make social decisions for the Herd
    def social_decision(self, world, mappa):
        if len(self.ErbastList) > 1:
            # Iterate over members starting from the second one
            for member in self.ErbastList[1:]:
                # Calculate the probability of removing the member based on its social aggression (sa)
                probability = 0.75 + (member.get_sa() / 20)
                if rdm.random() > probability:
                    # Remove the member and either add it to an existing group or create a new Herd
                    self.remove_and_create_herd(member, world[self.position][Erbast.__name__])
        # Move the Herd or stay in the current position
        self.move_or_stay(world, mappa)

    # Remove a member from the Herd and decide whether to add it to an existing group or create a new Herd
    def remove_and_create_herd(self, member, group_list):
        self.remove_member(member)
        if rdm.random() < member.get_sa() / 5:
            rdm.choice(group_list).add_erbast(member)
        else:
            group_list.append(Initialize_Herd(self.position, member))
        return

    # Move the Herd or stay in the current position based on strategy and world conditions
    def move_or_stay(self, world, mappa):
        avg_en = self.get_total_energy() / self.get_number_of_members()
        if self.strategy == 'Foraging':
            if world[self.position][Vegetob.__name__]:
                probability = 0.2 + (0.8 - (avg_en / 62.5))
            else:
                probability = 0
        elif self.strategy == 'Migration':
            if world[self.position][Vegetob.__name__]:
                probability = 1 - (avg_en / 50)
            else:
                probability = 0
        else:
            if world[self.position][Vegetob.__name__]:
                probability = 0.1 + (0.9 - (avg_en / 55.5))
            else:
                probability = 0

        if rdm.random() < probability:
            self.set_color("Green")  # Set color to "Green" and stay
        else:
            self.move_herd(world, mappa)  # Move the Herd
            self.set_color("Black")  # Set color to "Black" after movement
            self.set_color("Green")  # Set color to "Green" after movement

    # Allow the Herd to graze, removing it from the world if it has no members
    def graze(self, world):
        if not self.get_members_list():
            world[self.position][Erbast.__name__].remove(self)
            return
        for creature in self.get_members_list():
            if isinstance(creature, Erbast):
                creature.graze(world)
        return

    # Check the number of members in the Herd and remove excess members
    def check_members(self, group_list):
        while True:
            if self.get_number_of_members() > 15:
                self.remove_and_create_herd(rdm.choice(self.get_members_list()), group_list)
            else:
                break
        for member in self.get_members_list():
            member.check_conditions()
            if member.get_energy() <= 0:
                self.remove_member(member)
        return

    # Check if the Herd is empty (has no members)
    def check_if_empty(self):
        return not self.ErbastList

    # Allow the Herd to spawn new Erbast members
    def Spawn(self):
        for creature in self.get_members_list():
            if isinstance(creature, Erbast):
                creature.add_age()
                if creature.age_check():
                    self.remove_member(creature)
                    num_children = rdm.randint(1, 3)
                    for n in range(num_children):
                        self.add_erbast(
                            Initialize_Erbast(creature.position, None))
        return


# Define the Pride class
class Pride:
    def __init__(self, position):
        # Initialize a Pride with a target, strategy, position, and attributes
        self.target = None  # Initialize target as None
        self.strategy = rdm.choice(Carvix_Strategies + ['None'])  # Randomly select a strategy
        self.position = position  # Set the initial position
        self.CarvixList = []  # Initialize a list to store Carvix instances
        self.color = "White"  # Initialize the color attribute as "White"
        self.direction = None  # Initialize direction as None
        self.track = False  # Initialize track as False

    # Getter method to retrieve the direction
    def get_direction(self):
        return self.direction

    # Setter method to set the direction
    def set_direction(self, d):
        self.direction = d
        return

    # Check if a direction is set
    def check_direction(self):
        if self.direction:
            return True
        else:
            return False

    # Get the total energy of the Pride by summing the energy of its members
    def get_total_energy(self):
        return sum(member.get_energy() for member in self.get_members_list() or [])

    # Getter method to retrieve the color
    def get_color(self):
        return self.color

    # Setter method to set the color
    def set_color(self, color):
        self.color = color
        return

    # Reset the color to "White"
    def reset_color(self):
        self.color = "White"
        return

    # Add a Carvix instance to the Pride
    def add_carvix(self, carvix):
        if isinstance(carvix, Carvix):
            self.CarvixList.append(carvix)
        return

    # Get the number of Carvix members in the Pride
    def get_number_of_members(self):
        return len(self.CarvixList)

    # Get the list of Carvix members
    def get_members_list(self):
        return self.CarvixList

    # Remove a member from the Pride
    def remove_member(self, member):
        if member in self.CarvixList:
            self.CarvixList.remove(member)

    # Get the leader of the Pride based on energy and sa (social aggression)
    def get_leader(self):
        if self.CarvixList:
            lead = max(self.CarvixList, key=lambda x: x.energy + x.sa)
            self.CarvixList.remove(lead)
            self.CarvixList.insert(0, lead)
            return lead

    # Get the weakest member of the Pride based on energy
    def get_weaker(self):
        if self.CarvixList:
            weak = min(self.CarvixList, key=lambda x: x.energy)
            self.CarvixList.remove(weak)
            self.CarvixList.insert(-1, weak)
            return weak

    # Move the Pride in the world based on certain conditions
    def move_pride(self, world, mappa):
        # Store the old position
        old_position = self.position
        directions = [(1, 0), (-1, 0), (0, 1), (0, -1), (0, 0)]  # Possible directions
        weighted_directions = []  # Store valid directions

        # Iterate over possible directions
        for direction in directions:
            new_position = (old_position[0] + direction[0], old_position[1] + direction[1])

            # Check if the new position is within world boundaries
            if not (0 <= new_position[0] < world_size[0] and 0 <= new_position[1] < world_size[1]):
                continue

            # Check if the ground at the new position is suitable
            if Ground_Check(new_position, mappa):
                weighted_directions.append(direction)

        # If there are valid directions, select the one with the highest evaluation
        if weighted_directions:
            selected_direction = max(weighted_directions, key=lambda d: self.evaluate_direction(d))

            # If the selected direction is (0, 0), return
            if selected_direction == (0, 0):
                return

            # Update the position and move the Pride
            new_position = (old_position[0] + selected_direction[0], old_position[1] + selected_direction[1])
            if self in world[self.position][Carvix.__name__]:
                world[self.position][Carvix.__name__].remove(self)
            self.set_direction(selected_direction)
            self.position = new_position
            world[new_position][Carvix.__name__].append(self)

            # Update the position of the Carvix members in the Pride
            for c in self.CarvixList:
                c.rem_energy(MOVEMENT_COST)
                c.set_position(new_position)
        return

    # Evaluate the desirability of a direction for the Pride's movement
    def evaluate_direction(self, direction):
        weighted_sum = 0
        # Calculate the new position based on the direction
        new_position = (self.position[0] + direction[0], self.position[1] + direction[1])

        # Check if the new position is within the world boundaries
        if not (0 <= new_position[0] < world_size[0] and 0 <= new_position[1] < world_size[1]):
            return weighted_sum

        danger_level = 0
        # Iterate over a neighborhood around the new position
        for i in range(self.position[0] + direction[0] - Awareness, self.position[0] + direction[0] + Awareness + 1):
            for j in range(self.position[1] + direction[1] - Awareness,
                           self.position[1] + direction[1] + Awareness + 1):
                if not (0 <= i < world_size[0] and 0 <= j < world_size[1]):
                    continue
                if Ground_Check((i, j), mappa):
                    danger_level += calculate_dangerousness((i, j))
                    if not self.target:
                        self.target = (i, j)
                    if calculate_conditions_carvix((i, j)) > calculate_conditions_carvix((i, j)):
                        self.target = (i, j)

        energy_expenditure = calculate_energy_expenditure(calculate_distance(self.position, self.target), self)
        target_density = calculate_conditions_carvix(self.target) * 1000

        # Group Behavior: Factor in the proximity of other Pride members
        group_energy = calculate_energy_available(self)
        vegetob_density = calculate_conditions_carvix(new_position)

        if direction == self.direction:
            dir_weight = 1
        elif direction == opposite_direction(self.direction):
            dir_weight = 0.25
        else:
            dir_weight = 0.75

        if self.strategy == 'Fight':
            danger_level = danger_level * 0.5

        # Combine all factors with appropriate weights to compute a weighted sum
        weighted_sum = (
                dir_weight * target_density * vegetob_density * group_energy / energy_expenditure * danger_level
        )
        return weighted_sum

    # Make social decisions for the Pride members
    def social_decision(self, world, mappa):
        if not self.get_members_list():
            return
        if len(self.get_members_list()) == 1:
            self.move_or_stay(world, mappa)  # Make a decision for a single member Pride
            return
        for member in self.get_members_list():
            if self.strategy == 'Divide':
                probability = 0.35 + 0.05 * self.get_number_of_members()
            else:
                probability = 0.05 * self.get_number_of_members()
            if rdm.random() > probability:
                self.remove_and_create_pride(member,
                                             world[self.position][Carvix.__name__])  # Remove and create Pride members
        if self.get_members_list():
            self.move_or_stay(world, mappa)  # Move or stay after social decisions
        return

    # Remove a member and either add it to an existing Pride or create a new one
    def remove_and_create_pride(self, member, group_list):
        self.remove_member(member)
        if rdm.random() < 0.7:
            rdm.choice(group_list).add_carvix(member)
        else:
            group_list.append(Initialize_Pride(self.position, member))
        return

    # Move or stay based on the Pride's strategy and energy conditions
    def move_or_stay(self, world, mappa):
        if self.strategy == 'Insist':
            probability = 0.6 + (4 / self.get_total_energy())
        else:
            probability = 0.25 + (4 / self.get_total_energy())
        if rdm.random() < probability:
            self.set_color("Green")  # Set color to "Green" and stay
            return
        else:
            self.move_pride(world, mappa)  # Move the Pride
            self.set_color("Black")  # Set color to "Black" after movement
            return

    # Check the number of members in the Pride and remove excess members
    def check_members(self, group_list):
        while self.get_number_of_members() > 10:
            self.remove_and_create_pride(rdm.choice(self.get_members_list()), group_list)
        for member in self.get_members_list():
            if member.get_energy() <= 0:
                self.remove_member(member)
        return

    # Check if the Pride is empty (has no members)
    def check_if_empty(self):
        if not self.get_members_list():
            return True

    # Allow the Pride to spawn new Carvix members
    def Spawn(self):
        for creature in self.get_members_list():
            if isinstance(creature, Carvix):
                creature.add_age()
                if creature.age_check():
                    self.remove_member(creature)
                    num_children = rdm.randint(1, 3)
                    for n in range(num_children):
                        self.add_carvix(
                            Initialize_Carvix(creature.position, None))
        return

    # Hunt other herds in the world
    def hunt(self, herds_list):
        if self.check_if_empty():
            return
        target_herd = rdm.choice(herds_list)
        if not isinstance(target_herd, Herd):
            return
        while not target_herd.get_members_list():
            if not herds_list:
                return
            herds_list.remove(target_herd)
            target_herd = rdm.choice(herds_list)

        avg_energy_self = self.get_total_energy() / self.get_number_of_members()
        avg_energy_herd = target_herd.get_total_energy() / target_herd.get_number_of_members()

        # If the Pride has higher average energy, hunt the target herd
        if avg_energy_self > avg_energy_herd:
            while avg_energy_self > avg_energy_herd:
                c = rdm.choice(target_herd.get_members_list())
                if self.get_leader().get_energy() > c.get_energy():
                    self.get_leader().rem_energy(MOVEMENT_COST)
                    self.get_weaker().add_energy(HUNTING_ENERGY_GAIN)
                    for car in self.get_members_list():
                        car.add_energy(2)
                    target_herd.remove_member(c)
                else:
                    self.get_leader().rem_energy(MOVEMENT_COST * 3)
                if target_herd.check_if_empty() or target_herd.get_number_of_members() < 3:
                    return
                avg_energy_self = self.get_total_energy() / self.get_number_of_members()
                avg_energy_herd = target_herd.get_total_energy() / target_herd.get_number_of_members()

        # If the target herd has higher average energy, attempt to hunt a member
        else:
            c = rdm.choice(target_herd.get_members_list())
            if self.get_leader().get_energy() > c.get_energy():
                self.get_leader().rem_energy(MOVEMENT_COST)
                self.get_weaker().add_energy(HUNTING_ENERGY_GAIN)
                for car in self.get_members_list():
                    car.add_energy(2)
                target_herd.remove_member(c)
            else:
                self.get_leader().rem_energy(MOVEMENT_COST * 3)
        return


# Define the creature super class
class Creature:
    def __init__(self, position):
        self.position = position


# Define the Vegetob class, which inherits from the Creature class
class Vegetob(Creature):
    def __init__(self, position, density):
        # Initialize a Vegetob with a position, density, and attributes
        super().__init__(position)  # Call the constructor of the base class (Creature)
        self.density = density  # Set the initial density
        self.days_passed_since_growing = 0  # Initialize the number of days passed since growing

    # Method to simulate growth of Vegetob
    def grow(self):
        # Calculate the probability of growth based on the number of days passed since last growth
        probability = 0.3 + (self.get_days_passed_since_growing() * 0.07)

        # Check if the random number falls below the growth probability
        if rdm.random() < probability:
            # Increase the density by a random value between 3 and 5
            self.add_density(rdm.randint(3, 5))
        else:
            # Increment the number of days passed since growing
            self.days_passed_since_growing = self.get_days_passed_since_growing() + 1

    # Method to add density to the Vegetob
    def add_density(self, n):
        self.density = self.density + n
        return

    # Method to remove density from the Vegetob
    def remove_density(self, n):
        self.density = self.density - n
        return

    # Getter method to retrieve the density of the Vegetob
    def get_density(self):
        return self.density

    # Getter method to retrieve the number of days passed since growing
    def get_days_passed_since_growing(self):
        return self.days_passed_since_growing


# Define the Erbast class, which inherits from the Creature class
class Erbast(Creature):
    def __init__(self, position, energy, lifetime, age, sa):
        # Initialize an Erbast with a position, energy, lifetime, age, and sa (social ability)
        super().__init__(position)  # Call the constructor of the base class (Creature)
        self.energy = energy  # Set the initial energy level
        self.lifetime = lifetime  # Set the maximum lifetime of the Erbast
        self.age = age  # Set the initial age
        self.sa = sa  # Set the social ability

    # Method to check and adjust Erbast conditions
    def check_conditions(self):
        # Ensure that energy does not exceed the maximum energy limit
        if not self.energy < MAX_ENERGY_ERBAST:
            self.energy = MAX_ENERGY_ERBAST
        return

    # Method to add energy to the Erbast
    def add_energy(self, n):
        self.energy = self.energy + n
        return

    # Method to remove energy from the Erbast
    def rem_energy(self, n):
        self.energy = self.energy - n
        return

    # Getter method to retrieve the energy level of the Erbast
    def get_energy(self):
        return self.energy

    # Getter method to retrieve the social ability (sa) of the Erbast
    def get_sa(self):
        return self.sa

    # Method to set the position of the Erbast
    def set_position(self, new_position):
        self.position = new_position

    # Method for Erbast to graze on Vegetob and gain energy
    def graze(self, world):
        if world[self.position][Vegetob.__name__]:
            # Remove density from the Vegetob and add energy to the Erbast
            world[self.position][Vegetob.__name__][0].remove_density(5)
            self.add_energy(GRAZING_ENERGY_GAIN)

    # Method to increase the age of the Erbast
    def add_age(self):
        self.age += 1
        return

    # Method to check if the Erbast has reached its maximum lifetime
    def age_check(self):
        if self.age == self.lifetime:
            return True
        else:
            return False


# Define the Carvix class, which inherits from the Creature class
class Carvix(Creature):
    def __init__(self, position, energy, lifetime, age, sa):
        # Initialize a Carvix with a position, energy, lifetime, age, and sa (social ability)
        super().__init__(position)  # Call the constructor of the base class (Creature)
        self.energy = energy  # Set the initial energy level
        self.lifetime = lifetime  # Set the maximum lifetime of the Carvix
        self.age = age  # Set the initial age
        self.sa = sa  # Set the social ability

    # Method to check and adjust Carvix conditions
    def check_conditions(self):
        # Ensure that energy does not exceed the maximum energy limit
        if self.energy > MAX_ENERGY_CARVIX:
            self.energy = MAX_ENERGY_CARVIX
        return

    # Method to set the energy level of the Carvix
    def set_energy(self, n):
        self.energy = n
        return

    # Method to add energy to the Carvix
    def add_energy(self, n):
        self.energy = self.energy + n
        return

    # Getter method to retrieve the social ability (sa) of the Carvix
    def get_sa(self):
        return self.sa

    # Method to remove energy from the Carvix
    def rem_energy(self, n):
        self.energy = self.energy - n
        return

    # Getter method to retrieve the energy level of the Carvix
    def get_energy(self):
        return self.energy

    # Method to set the position of the Carvix
    def set_position(self, new_position):
        self.position = new_position

    # Method to increase the age of the Carvix
    def add_age(self):
        self.age += 1
        return

    # Method to check if the Carvix has reached its maximum lifetime
    def age_check(self):
        if self.age == self.lifetime:
            return True
        else:
            return False


# Planissus Constants:

# Dimension of the world simulated
world_size = (80, 80)

# On every update, the number of new Vegetob that appear
TotGrowth = world_size[0]

# Number of cells in the world
total_cells = world_size[0] * world_size[1]

# Number of creatures that are present at the start if the simulation
num_creatures = world_size[0] ** 2

# Types of creatures
creature_types = [Vegetob, Erbast, Carvix]

# Number of cells that a creature can analyze
Awareness = 5

# Constants for maximum energy and density
MAX_ENERGY_ERBAST = 50
MAX_ENERGY_CARVIX = 60
MAX_DENSITY_VEGETOB = 100

# Constants for maximum lifetime
MAX_LIFETIME_ERBAST = 100
MAX_LIFETIME_CARVIX = 75

# Energy and movement-related constants
MOVEMENT_COST = 1
GRAZING_ENERGY_GAIN = 15
HUNTING_ENERGY_GAIN = 10

# Strategies for Erbast and Carvix
Erbast_Strategies = ['Migration', 'Foraging', 'Escape', 'Fight']
Carvix_Strategies = ['Insist', 'Divide', 'Fight']

# Colors for visualization
colors = ['#0278c7', '#589130', '#006400', '#FFB90F', '#cc2702', '#000B55']
cmap = ListedColormap(colors)

# Variable to track days passed
Days_Passed = 0

# Start of the Planissus program

world = World_Generation()
mappa = Map_generator()
World_Initialization()

# Create the main window
window = ttk.Window(themename='solar')
window.title("Planissus")

# Disable window resizing
window.resizable(False, False)

# Create a Frame to hold the entire content
main_frame = tk.Frame(window)
main_frame.pack(expand=True)  # Expand to fill the window

# Create a Frame for the plot
plot_frame = tk.Frame(main_frame, bg="gray14")
plot_frame.pack(side="left", padx=20, expand=True)  # Expand to fill the left side

# Create a Figure for the plot
fig = Figure(figsize=(6, 6))
fig.set_facecolor('#0278c7')
ax = fig.add_subplot()
ax.axis('off')
fig.canvas.mpl_connect('button_press_event', inspectCell)

# Add a title label
Guide_Label = tk.Label(plot_frame, text="Welcome to Planissus!", font=("Fixedsys", 20), fg="white", bg="gray14")

Guide_Label.pack(padx=100, pady=20, anchor="center")

# Create a canvas for the plot
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack(expand=True, fill="both")  # Expand to fill the plot_frame

mat = ax.matshow(Generate_Grid(world, mappa), cmap=cmap, vmin=0, vmax=5)
main_animation = FuncAnimation(fig, update, frames=100, interval=100)
anim_running = True

# Bind mouse motion events to update the cursor coordinates label
canvas.mpl_connect('motion_notify_event', update_cursor_label)

# Create a Frame for the labels with fixed dimensions
label_frame = tk.Frame(main_frame)  # Adjust width and height as needed
label_frame.pack(side="bottom")

# Create a label to display cursor coordinates
cursor_label = tk.Label(label_frame, text="", fg="white", bg="gray14", anchor="sw", width=20, height=2, wraplength=140)
cursor_label.grid(row=1, column=0, sticky="sw")

# Create a label to indicate days passed on the bottom right
days_passed_label = tk.Label(label_frame, text="Days Passed: 0", fg="white", bg="gray14", anchor="se")
days_passed_label.grid(row=1, column=1, sticky="se")

# Create a Frame for buttons with the desired background color
buttons_frame = tk.Frame(main_frame, bg="gray14")
buttons_frame.pack(side="right", padx=20)

# Style the buttons
Pause_Button = tk.Button(buttons_frame, text="Pause", command=press_pause_button, bg="gray", fg="white", padx=10,
                         pady=5, font=('Fixedsys', 20))
Restart_Button = tk.Button(buttons_frame, text="Restart Simulation", command=restart_sim_button, bg="gray", fg="white",
                           padx=10, pady=5, state="disabled", font=('Fixedsys', 15))
New_Map_Button = tk.Button(buttons_frame, text="New Map", command=new_map_button, bg="gray", fg="white", padx=10,
                           pady=5, state="disabled", font=('Fixedsys', 15))

# Pack the buttons with increased spacing
Pause_Button.pack(pady=10)
Restart_Button.pack(pady=10)
New_Map_Button.pack(pady=10)

# Start the GUI main loop
window.mainloop()