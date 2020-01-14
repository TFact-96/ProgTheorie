from test2 import Lattice2DEnv
from gym import spaces
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import pandas as pd

LEFT_CMD = [1, 0, 0, 0]
DOWN_CMD = [0, 1, 0, 0]
UP_CMD = [0, 0, 1, 0]
RIGHT_CMD = [0, 0, 0, 1]

MIN_REWARD = 70

np.random.seed(42)

seq = 'HHPHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'  # Our input sequence
action_space = spaces.Discrete(4)  # Choose among [0, 1, 2 ,3]
env = Lattice2DEnv(seq)


def play_random_games(games=100):
    """
    Play Random Games to Get Some Observations
    :param games:
    :return:
    """

    # Storage for All Games Movements
    all_movements = []

    for episode in range(games):

        # Reset Game Reward
        episode_reward = 0

        # Define Storage for Current Game Data
        current_game_data = []

        # Reset Game Environment
        env.reset()

        # Get First Random Movement
        action = action_space.sample()

        while True:
            # Play
            observation, reward, done, info = env.step(action)

            # Get Random Action (On Real, its get a "Next" movement to compensate Previous Movement)
            action = env.action_space.sample()

            if action == 0:
                direction = LEFT_CMD
            elif action == 1:
                direction = DOWN_CMD
            elif action == 2:
                direction = UP_CMD
            else:
                direction = RIGHT_CMD

            # Store Observation Data and Action Taken
            current_game_data.append(
                np.hstack((observation, direction))
            )

            if done:
                break

            episode_reward -= reward

        # Save All Data (Only for the Best Games)
        if episode_reward >= MIN_REWARD:
            print('.', end='')
            all_movements.extend(current_game_data)

    # Create DataFrame
    dataframe = pd.DataFrame(
        all_movements,
        columns=['chain_length', 'seq_length', 'collisions', 'is_trapped', 'action_to_left', 'action_to_down',
                 'action_to_up', 'action_to_right']
    )
    print(dataframe)

    # Convert Action Columns to Integer
    dataframe['action_to_left'] = dataframe['action_to_left'].astype(int)
    dataframe['action_to_down'] = dataframe['action_to_down'].astype(int)
    dataframe['action_to_up'] = dataframe['action_to_up'].astype(int)
    dataframe['action_to_right'] = dataframe['action_to_right'].astype(int)

    return dataframe


def generate_ml(dataframe):
    """
    Train and Generate NN Model
    :param dataframe:
    :return:
    """

    # Define Neural Network Topology
    model = Sequential()
    model.add(Dense(64, input_dim=4, activation='relu'))
    # model.add(Dense(128,  activation='relu'))
    # model.add(Dense(128,  activation='relu'))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(4, activation='sigmoid'))

    # Compile Neural Network
    model.compile(optimizer='adam', loss='categorical_crossentropy')

    # Fit Model with Data
    model.fit(
        dataframe[['chain_length', 'seq_length', 'collisions', 'is_trapped']],
        dataframe[['action_to_left', 'action_to_down', 'action_to_up', 'action_to_right']],
        epochs=20
    )

    return model


def play_game(ml_model, games=100):
    """
    Play te Game
    :param ml_model:
    :param games:
    :return:
    """

    for i_episode in range(games):

        # Define Reward Var
        episode_reward = 0

        # Reset Env for the Game
        observation = env.reset()

        while True:
            #env.render()
            # Predict Next Movement
            current_action_pred = ml_model.predict(observation.reshape(1, 4))

            # Define Movement
            current_action = np.argmax(current_action_pred)
            print(current_action)

            # Make Movement
            observation, reward, done, info = env.step(current_action)

            if done:
                episode_reward += 1
                # print(f"Episode finished after {i_episode + 1} steps", end='')
                break

            episode_reward += 1

        print(f" Score = {episode_reward}")


print("[+] Playing random games")
df = play_random_games(games=1000)

print("[+] Training NN Model")
ml_model = generate_ml(df)

print("[+] Playing Games with NN")
play_game(ml_model=ml_model, games=100)
