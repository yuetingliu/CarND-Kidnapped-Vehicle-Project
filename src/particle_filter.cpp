/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#include "map.h"

using std::string;
using std::vector;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  // Set the number of particles
  num_particles = 1000;
  weights.assign(num_particles, 1);

  // std::vector<Particle> particles (num_particles);

  // random engine
  std::default_random_engine gen;
  // set standard deviations
  double std_x, std_y, std_theta;
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];

  // create a normal distribution for x, y, and theta
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);

  double sample_x, sample_y, sample_theta;
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;
  double weight = 1;
  for (int i=0; i<num_particles; ++i) {

      sample_x = dist_x(gen);
      sample_y = dist_y(gen);
      sample_theta = dist_theta(gen);

      Particle p{
        i, sample_x, sample_y,
        sample_theta, weight,
        associations, sense_x, sense_y
      };
      particles.push_back(p);

  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  for (int i=0; i<num_particles; ++i) {
      Particle p = particles[i];
      // update x
      double x_change = sin(p.theta+yaw_rate*delta_t) - sin(p.theta);
      double x_new = p.x + velocity / yaw_rate * x_change;
      // update y
      double y_change = cos(p.theta) - cos(p.theta+yaw_rate*delta_t);
      double y_new = p.y + velocity / yaw_rate * y_change;
      // update theta
      double theta_new = p.theta + yaw_rate * delta_t;

      // add Gaussian noise
      std::default_random_engine gen;
      // set standard deviations
//      double std_x, std_y, std_theta;
//      std_x = std_pos[0];
//      std_y = std_pos[1];
//      std_theta = std_pos[2];

//      // create a normal distribution for x, y, and theta
//      std::normal_distribution<double> dist_x(x_new, std_x);
//      std::normal_distribution<double> dist_y(y_new, std_y);
//      std::normal_distribution<double> dist_theta(theta_new, std_theta);

      // update particle
//      p.x = dist_x(gen);
//      p.y = dist_y(gen);
//      p.theta = dist_theta(gen);
//      particles[i] = p;
      particles[i].x = x_new;
      particles[i].y = y_new;
      particles[i].theta = theta_new;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    for (int i=0; i<num_particles; i++) {
        Particle p = particles[i];
        double x_part = p.x;
        double y_part = p.y;
        double theta_part = p.theta;
        std::vector<int> p_associations;
        std::vector<double> p_sense_x, p_sense_y;

        double new_weight = 1;
        for (unsigned int j=0; j<observations.size(); j++) {
            LandmarkObs land_obs = observations[j];
            double x_obs = land_obs.x;
            double y_obs = land_obs.y;

            // transform observation from vehicle coordinate to map coordinate
            double x_map;
            x_map = x_part + (cos(theta_part) * x_obs) - (sin(theta_part) * y_obs);
            double y_map;
            y_map = y_part + (sin(theta_part) * x_obs) + (cos(theta_part) * y_obs);

            // associate observation landmarks to map landmarks
            int size_landmark = map_landmarks.landmark_list.size();
            double diff;
            double min_diff = sensor_range;  // initialize min diff
            Map::single_landmark_s asso_landmark, current_landmark;
            for (int k=0; k<size_landmark; k++) {
                current_landmark = map_landmarks.landmark_list[k];
                float x_landmark = current_landmark.x_f;
                float y_landmark = current_landmark.y_f;
                diff = dist(x_landmark, y_landmark, x_map, y_map);
                if (diff < min_diff) {
                    min_diff = diff;
                    asso_landmark = current_landmark;
                }
            }

            //update particle association
            int id_i = asso_landmark.id_i;
            p_associations.push_back(id_i);
            p_sense_x.push_back(x_map);
            p_sense_y.push_back(y_map);

            // get associated landmark coordinates
            double mu_x = asso_landmark.x_f;
            double mu_y = asso_landmark.y_f;
            double std_x = std_landmark[0];
            double std_y = std_landmark[1];

            // calculate weights and update particle weight
            double weight = multivar_prob(std_x, std_y, x_map, y_map,
                                          mu_x, mu_y);

            new_weight *= weight;
        }
        particles[i].weight = new_weight;
        SetAssociations(particles[i], p_associations, p_sense_x, p_sense_y);
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

    std::default_random_engine generator;
    // update weights
    for (int i; i<num_particles; i++) {
        weights[i] = particles[i].weight;
    }

    std::discrete_distribution <int> distribution (weights.begin(), weights.end());

    std::vector <Particle> new_particles (num_particles);
    for (int i=0; i<num_particles; i++) {
       int p_ind = distribution(generator);
       new_particles[i] = particles[p_ind];
    }
    particles = new_particles;
}

double ParticleFilter::multivar_prob(double std_x, double std_y, double x_obs,
                                     double y_obs, double mu_x, double mu_y) {
    // first part of gaussain
    double gauss_part1;
    gauss_part1 = 1 / (2 * M_PI * std_x * std_y);

    // exponential part of gaussain
    double exponent;
    exponent = (pow(x_obs - mu_x, 2) / (2 * pow(std_x, 2))) +
               (pow(y_obs - mu_y, 2) / (2 * pow(std_y, 2)));

    // get weigth
    double weight;
    weight = gauss_part1 * exp(-exponent);

    return weight;
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
