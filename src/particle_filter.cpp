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
#include <limits.h>
#include  <cstdint>
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
// NOTE: Consult particle_filter.h for more information about this method (and others in this file).:q

	default_random_engine gen;
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	num_particles = 10000;

	particles = vector<Particle>(num_particles);

	for (int i = 0; i < num_particles; ++i) {
		Particle pt;
		pt.id = i;
		pt.x = dist_x(gen);
		pt.y = dist_y(gen);
		pt.theta = dist_theta(gen);
		pt.weight = 0.0;
		particles[i] = pt;
	}


	weights = vector<double>(num_particles, 0.0);
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);


	for (int i = 0; i < particles.size(); ++i) {
		particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(gen);
		particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) + dist_y(gen);
		particles[i].theta += particles[i].theta + yaw_rate * delta_t + dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < observations.size(); ++i) {
		double min_dist = INT_MAX;
		int min_id = -1;
		for (int j = 0; j < predicted.size(); ++j) {
			double curr_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (curr_dist < min_dist) {
				min_dist = curr_dist;
				min_id = j;
			}
		}

		observations[i].id = min_id;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	vector<LandmarkObs> land_marks = vector<LandmarkObs>(map_landmarks.landmark_list.size());
	for (int i = 0; i < land_marks.size(); ++i) {
		land_marks[i].id = map_landmarks.landmark_list[i].id_i;
		land_marks[i].x = map_landmarks.landmark_list[i].x_f;
		land_marks[i].y = map_landmarks.landmark_list[i].y_f;
	}

	for (int i = 0; i < particles.size(); ++i) {
		std::vector<LandmarkObs> transformed_obs = observations;
		double  weight = 1.;
		for (int j = 0; j < observations.size(); ++j) {
			// Transform car's observations to map coordinate system
			transformed_obs[j].x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
			transformed_obs[j].y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
		}
		
		dataAssociation(transformed_obs, land_marks);
		
		for (int j = 0; j < land_marks.size(); ++j) {
			int id = land_marks[j].id;
			double landmark_dist = dist(particles[i].x, particles[i].y, land_marks[j].x, land_marks[j].y);
			if (landmark_dist <=  sensor_range){
				double xmu = transformed_obs[id].x - land_marks[j].x;
				double ymu = transformed_obs[id].y - land_marks[j].y;
				double sigmax = std_landmark[0];
				double sigmay = std_landmark[1];				
				double w = exp(-((xmu * xmu / (2 * sigmax*sigmax)) + (ymu * ymu / (2 * sigmay*sigmay)))) / (2 * M_PI * sigmax * sigmay);
				weight *= w;				
			}
		}

		particles[i].weight = weight + 0.0001;
		weights[i] = weight + 0.0001;
	    
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> new_particles = vector<Particle>(particles.size());
	random_device rd;
	mt19937 gen(rd());
	discrete_distribution<> d (weights.begin(), weights.end());
	for (int n = 0; n < particles.size(); ++n) {
		 int idx = d(gen);
		 new_particles[n] = particles[idx];
	}
	particles = new_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
