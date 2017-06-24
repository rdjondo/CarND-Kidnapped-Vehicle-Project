/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};



class ParticleFilter {
	
	// Number of particles to draw
	int num_particles; 

	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	std::vector<double> weights;



public:
	
	// Set of current particles
	std::vector<Particle> particles;

	// Constructor
	// @param M Number of particles
	ParticleFilter() : num_particles(1), is_initialized(false) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);


	/**
	 * Compute real world position coordinates of the landmarks from the point of view of each particle
	 * @param observations
	 * @param sensor_range
	 * @param map_landmarks
	 */
	void positionsToWorld(std::vector<LandmarkObs> &observations, Map &map_landmarks);

    /**
     * Calculate position of average Particle
     * @return
     */
    Particle calculateMeanParticle() const;

    /**
     * Gather the world coordinates landmark position inside the 6-sigma circle plus
     * sensor range around average particle
     * @param sensor_range
     * @param map_landmarks
     * @param mean_particle
     * @param predicted
     */
    void findCandidateLandMarks(double sensor_range, const Map &map_landmarks, const Particle &mean_particle,
                                std::vector<LandmarkObs> &predicted) const;

	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations,
						 double sensor_range, Map &map_landmarks);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
	 *   standard deviation of bearing [rad]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations,
			Map map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();


	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}

    void nearestNeighbour(std::vector<LandmarkObs> &predicted, Map &map_landmarks);
};



#endif /* PARTICLE_FILTER_H_ */
