/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles.
	num_particles = 3;
	particles.reserve(num_particles);
	weights.reserve(num_particles);

	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	//Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// We will use 2 independent Gaussian distributions in first approximation
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; ++i) {
        Particle particle = {};
		particle.x = dist_x(gen);
		particle.id = i;
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1;
        particles.push_back(particle);
	}


	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
		double velocity, double yaw_rate) {
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];

	normal_distribution<double> noise_x(0.0, std_x);
	normal_distribution<double> noise_y(0.0, std_y);
	normal_distribution<double> noise_theta(0.0, std_theta);

	// Add measurements to each particle to integrate the position using velocity and yaw rate
	for (int i = 0; i < num_particles; ++i) {
		double xf;
		double yf;
		double theta_f;

		if(fabs(yaw_rate)<1e-4){
			xf = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			yf = particles[i].x + velocity * delta_t * sin(particles[i].theta);
			theta_f = particles[i].theta;
		} else {
			theta_f = fmod(particles[i].theta + yaw_rate * delta_t, 2*M_PI);
			xf = particles[i].x + velocity / yaw_rate * (sin(theta_f) - sin(particles[i].theta));
			yf = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(theta_f));
		}

		// Add random Gaussian noise
		particles[i].x = xf + noise_x(gen);
		particles[i].y = yf + noise_y(gen);
		particles[i].theta = fmod(theta_f + noise_theta(gen), 2*M_PI);
	}
}


void ParticleFilter::positionsToWorld(std::vector<LandmarkObs> &observations, Map &map_landmarks){

	for(int i = 0; i < num_particles; ++i){
		cout<<" particle id:"<<particles[i].id<<"  x:"<< particles[i].x << "  y:" << particles[i].y <<endl;

		particles[i].associations.clear();
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();

		particles[i].associations.reserve(observations.size());
		particles[i].sense_x.reserve(observations.size());
		particles[i].sense_y.reserve(observations.size());

		double theta_to_world = -particles[i].theta;
		cout<<"Transform observations to world coordinates :"<<endl;

		// Transform sensed landmarks in world coordinates as if the data was sensed from each particle.
		for (int j = 0; j < observations.size(); ++j) {
			particles[i].sense_x.push_back( observations[j].x * cosf(theta_to_world)
											- observations[j].y * sinf(theta_to_world) - particles[i].x );
			particles[i].sense_y.push_back( observations[j].x * sinf(theta_to_world)
											+ observations[j].y * cosf(theta_to_world) - particles[i].y );
			cout << "obs id:" << j << "  x:"
				 << particles[i].sense_x[j] << "  y:" << particles[i].sense_y[j]
				 << endl;
		}

		cout<<"Particle sensed size vector :"<< particles[i].sense_x.size()<<endl;
	}
}

Particle ParticleFilter::calculateMeanParticle() const {
	Particle mean_particle;
	mean_particle.x=0.0;
	mean_particle.y=0.0;
	for(int i = 0; i < num_particles; ++i){
		mean_particle.x += particles[i].x;
		mean_particle.y += particles[i].y;
	}
	if(num_particles>0) {
		mean_particle.x = mean_particle.x / num_particles;
		mean_particle.y = mean_particle.y / num_particles;
	} else{
		cout<<"ERROR: Number of particles should be greater than 0"<<endl;
	}

    cout<<" Mean particle x:"<< mean_particle.x << "  y:" << mean_particle.y <<endl;
	return mean_particle;
}


void ParticleFilter::findCandidateLandMarks(double sensor_range, const Map &map_landmarks, const Particle &mean_particle,
											vector<Map::single_landmark_s> &candidate_landmarks) const {
	const double six_sigma = 6 * sqrt(std_x * std_x + std_y * std_y);
	for (Map::single_landmark_s landmark: map_landmarks.landmark_list) {
		if (dist(landmark.x_f, landmark.y_f, mean_particle.x, mean_particle.y) <
			sensor_range + six_sigma) {
			Map::single_landmark_s landmark_candidate = landmark;
			candidate_landmarks.push_back(landmark_candidate);
		}
	}

	cout<<"candidate_landmarks.size()="<<candidate_landmarks.size()<<endl;
	if(candidate_landmarks.size()==0){
		cout<<"ERROR: Could not find landmarks around particles"<<endl;
	}
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations,
                                     double sensor_range, Map &map_landmarks) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    cout<<" after positionsToWorld Count all particles , number : "<<particles.size()<<endl;

	// Compute real world position coordinates of the landmarks from the point of view of each particle
	positionsToWorld(observations, map_landmarks);

    cout<<" after positionsToWorld Count all particles , number : "<<particles.size()<<endl;

    // Get mean position of all particles
	Particle mean_particle = calculateMeanParticle();

    cout<<" after calculateMeanParticle Count all particles , number : "<<particles.size()<<endl;

	// List all candidate landmarks around sensor range + 6-sigma of mean particle
	// The idea is to reduce the number of landmarks to compare to the sensed landmarks
    predicted.clear();
    predicted.reserve(5);

    std::vector<Map::single_landmark_s> candidate_landmarks;
    candidate_landmarks.clear();
    candidate_landmarks.reserve(11);
	findCandidateLandMarks(sensor_range, map_landmarks, mean_particle, candidate_landmarks);

    cout<<" asgff Count all particles , number : "<<particles.size()<<endl;
	for(Particle particle : particles){
		// For each particle
        cout<<" Loop accross particles "<<endl;
        particle.associations.reserve(particle.sense_x.size());
        particle.associations.clear();
        std::vector<Map::single_landmark_s> candidate_landmarks_buf(candidate_landmarks);

        for (int i = 0; i < particle.sense_x.size(); ++i) {
			// compute the distance between each landmark position and the particule's sensed positions
            cout<<" Loop accross sensed position "<<endl;
			double min_dist_landmark = 1e6;
            int j;
            for ( j = 0; j < candidate_landmarks_buf.size(); ++j) {
				// Retain the association for the closest landmark
                auto candidate_landmark = candidate_landmarks_buf[j];
				double distance = dist(candidate_landmark.x_f, candidate_landmark.y_f, particle.sense_x[i],
									   particle.sense_y[i]);
				if (distance < min_dist_landmark) {
					min_dist_landmark = distance;
					particle.associations[i] = candidate_landmark.id_i;
				}
			}
            // Optimization :
            candidate_landmarks_buf.erase(candidate_landmarks_buf.begin()+j-1);

            cout<<"Particle "<<particle.id<<" sensed landmark "<<i<< " at distance "<< min_dist_landmark << " associated with landmark "<<particle.associations[i]<<endl;
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    std::vector<LandmarkObs> predicted;
    dataAssociation(predicted, observations, sensor_range, map_landmarks);

    // Get particle-measurements probability as a function of the distance to its associated landmark
    auto landmark_list = map_landmarks.landmark_list;
    const double sigma = sqrt(std_x * std_x + std_y * std_y);
    for(Particle particle:particles){
        double prob = 1.0;
        for (int i = 0; i < particle.associations.size(); ++i) {
            int associated_landmark = particle.associations[i];
            Map::single_landmark_s landmark = landmark_list[associated_landmark];
            double distance_predicted = dist(particle.x, particle.y,  landmark.x_f, landmark.y_f);
            double distance_measured = dist(observations[i].x, 0, observations[i].y, 0);
            prob *= unimodal_gaussian(distance_predicted, sigma, distance_measured);
        }
        particle.weight = prob;
        cout<< "particle.weight :"<<particle.weight<<endl;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle,
		std::vector<int> associations, std::vector<double> sense_x,
		std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
