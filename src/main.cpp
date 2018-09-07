#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// update car speed
void updateVelocity(double &vel, double tgt_vel, int state, double offset, bool car_close_to_front)
{
	double max_vel = 49.5;
	
	if(state==2)
	{
		if(vel < tgt_vel)
		{
			vel += .224;
		}
			
		if(!car_close_to_front)
		{
			tgt_vel = max_vel;
		}
	}
	else
	{
		if(car_close_to_front)
		{
			if(vel > tgt_vel)
			{
				vel -= (.224 + offset*1);
			}
			// else
			// {
				// tgt_vel -= offset*2;
			// }
		}
		else if(vel < 40.0)
		{
			vel += .224;
		}
	}
	
	// clamp reference velocity to max
	vel = (vel>max_vel)? max_vel: vel;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  
  int curr_lane = 1;			// Current car lane
  int tgt_lane = 1;				// Target car lane
  double curr_ref_vel = 0.0;	// Current reference velocity
  double tgt_ref_vel = 49.5; 	// Target reference velocity (mph)
  double close_dist_offset = 0.0;

  // 0 - maintaining current lane
  // 1 - preparing for lane change
  // 2 - cruising
  // 3 - making a lane change
  int car_state = 2;

  h.onMessage([
	&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
	&curr_lane,&curr_ref_vel,&tgt_lane,&tgt_ref_vel,&car_state,&close_dist_offset
  ](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
			
			int prev_size = previous_path_x.size();

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
			if(prev_size > 0)
			{
				car_s = end_path_s;
			}
			
			bool too_close_front = false;
			bool too_close_left = false;
			bool too_close_right = false;
			
			bool empty_lane_mid = true;
			bool empty_lane_left = true;
			bool empty_lane_right = true;
			
			int fastest_lane = 0;
			double fastest_speed = 0.0;
						
			// find ref_vel to use
			for(int i=0; i<sensor_fusion.size(); ++i)
			{
				float d = sensor_fusion[i][6];
				int lane = -1;
				
				double check_car_s = sensor_fusion[i][5];
				
				// // car detected by sensor is outside the range of concern
				// if(d<0 || ((check_car_s-car_s)<-10) || ((check_car_s-car_s)>30))
				// {
					// continue;
				// }
				
				double vx = sensor_fusion[i][3];
				double vy = sensor_fusion[i][4];
				double check_speed = sqrt(vx*vx + vy*vy);
				check_car_s += ((double)prev_size*.02*check_speed);
				
				// determine traffic conditions i.e. empty lanes and the fastest lane
				for(int i=0; i<3; ++i)
				{
					lane = i;
					if(d>(2+4*lane-2) && d<(2+4*lane+2) && 
					   (check_car_s-car_s)>0 && (check_car_s-car_s)<40)
					{
						if(i == 0)
							empty_lane_left = false;
						if(i == 1)
							empty_lane_mid = false;
						if(i == 2)
							empty_lane_right = false;
						
						if(check_speed > fastest_speed)
						{
							fastest_lane = i;
							fastest_speed = check_speed;
						}
					}
				}
				
				// determine if there are other cars close to my car
				for(int i=0; i<3; ++i)
				{
					lane = curr_lane + i-1;
					
					if(lane<0 || lane>2)
						continue;
					
					if(d>(2+4*lane-2) && d<(2+4*lane+2))
					{
						double left_of_lane = 0.0 - (12.0*abs(i-1));
						double right_of_lane = 30.0 - (20.0*abs(i-1));
						if(((check_car_s-car_s)>left_of_lane) && ((check_car_s-car_s)<right_of_lane))
						{
							if(i == 0)
								too_close_left = true;
							if(i == 1)
							{
								too_close_front = true;
								
								// adjust target reference velocity based on the distance and
								// the velocity of the car in front of my car
								close_dist_offset = check_car_s-car_s;
								close_dist_offset = (close_dist_offset>0)? 1/close_dist_offset: 0;
								tgt_ref_vel -= close_dist_offset*2;
								
								double speed = 20.0 + check_speed;
								tgt_ref_vel = (tgt_ref_vel<speed)? speed: tgt_ref_vel;
							}
							if(i == 2)
								too_close_right = true;
						}
					}
				}
			}
			
			// unset target lane if the car is currently cruising
			if(car_state == 2)
			{
				tgt_lane = -1;
				close_dist_offset = 0;
			}
			// do not update the target lane while the car is making a lane change
			else if(car_state == 3)
			{
				// do nothing
			}
			// if there is an empty lane available, attempt to update the target lane
			else if(empty_lane_left || empty_lane_mid || empty_lane_right)
			{
				// if the current lane is empty, no need to switch to the another empty lane
				if((curr_lane==0 && empty_lane_left) ||
				   (curr_lane==1 && empty_lane_mid) ||
				   (curr_lane==2 && empty_lane_right))
				{
					tgt_lane = -1;
				}
				// if current lane is not empty
				else
				{
					// TODO: change this to favor middle lane to reduce jerk motion across extremes
					if(empty_lane_left)
						tgt_lane = 0;
					else if(empty_lane_mid)
						tgt_lane = 1;
					else if(empty_lane_right)
						tgt_lane = 2;
				}
			}
			// if there is no empty lane available, unset the target lane
			else
			{
				tgt_lane = -1;
			}

			// debug prints
			if(false)
			{
				std::cout << "car_state: " << car_state << "    too_close: " << too_close_left 
						  << ", " << too_close_front << ", " << too_close_right 
						  << "    empty_lane: " << empty_lane_left << ", " << empty_lane_mid 
						  << ", " << empty_lane_right << "    lane_vars: " << curr_lane 
						  << ", " << tgt_lane << ", " << fastest_lane 
						  // << "    tgt_ref_vel: " << tgt_ref_vel
						  // << "    close_dist_offset: " << close_dist_offset
						  << std::endl;
			}

			updateVelocity(curr_ref_vel, tgt_ref_vel, car_state, close_dist_offset, too_close_front);
			
			// implement state machine for car's state
			if(car_state == 2)
			{	
				// if there is a car in my front
				// change my car's state from cruising to maintaining/following
				if(too_close_front)
				{
					car_state = 0;
				}
			}
			else if(car_state == 0)
			{
				// if there is a target lane and it is not equal to the current lane
				// change the car's state to prepare to make a lane change
				if((tgt_lane >= 0) && (tgt_lane != curr_lane))
				{
					car_state = 1;
				}
				else if(!too_close_front)
				{
					car_state = 2;
				}
			}
			else if(car_state == 1)
			{
				// if the target lane is invalid,
				// transition the car's state to maintaining the current lane
				if(tgt_lane < 0)
				{
					car_state = 0;
				}
				else
				{
					// wait until there are no other cars close to my car before
					// attempting to make lane changes 
					if(curr_lane > tgt_lane)
					{
						if(!too_close_left)
						{
							curr_lane--;
							car_state = 3;
						}
					}
					else if(curr_lane < tgt_lane)
					{
						if(!too_close_right)
						{
							curr_lane++;
							car_state = 3;
						}
					}
				}
			}
			else if(car_state == 3)
			{
				// if the target lane is invalid,
				// transition the car's state to maintaining the current lane
				if(tgt_lane < 0)
				{
					car_state = 0;
				}
				else
				{
					// wait until my car is within the center threshold of the current lane
					// before making other state transitions
					double threshold = 0.5;
					if(car_d>(2+4*curr_lane-threshold) && car_d<(2+4*curr_lane+threshold))
					{
						if(curr_lane != tgt_lane)
							car_state = 1;
						else
							car_state = 2;
					}
				}
			}
			
          	json msgJson;

          	vector<double> pts_x;
          	vector<double> pts_y;
			
			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);
			
			// if previous size is almost empty, use the car as starting reference
			if(prev_size < 2)
			{
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);
				
				pts_x.push_back(prev_car_x);
				pts_x.push_back(car_x);
				
				pts_y.push_back(prev_car_y);
				pts_y.push_back(car_y);
			}
			// use the previous path's and point as starting reference
			else
			{
				ref_x = previous_path_x[prev_size-1];
				ref_y = previous_path_y[prev_size-1];
				
				double ref_x_prev = previous_path_x[prev_size-2];
				double ref_y_prev = previous_path_y[prev_size-2];
				ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);
				
				pts_x.push_back(ref_x_prev);
				pts_x.push_back(ref_x);
				
				pts_y.push_back(ref_y_prev);
				pts_y.push_back(ref_y);
			}

			// In Frenet add evenly 30m spaced points ahead of the starting reference
			int spacing = 30;
			vector<double> next_wp0 = getXY(car_s+(spacing*1), (2+4*curr_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp1 = getXY(car_s+(spacing*2), (2+4*curr_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp2 = getXY(car_s+(spacing*3), (2+4*curr_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			
			pts_x.push_back(next_wp0[0]);
			pts_x.push_back(next_wp1[0]);
			pts_x.push_back(next_wp2[0]);
			
			pts_y.push_back(next_wp0[1]);
			pts_y.push_back(next_wp1[1]);
			pts_y.push_back(next_wp2[1]);
			
			for(int i=0; i<pts_x.size(); ++i)
			{
				// shift car reference angle to 0 degrees
				double shift_x = pts_x[i]-ref_x;
				double shift_y = pts_y[i]-ref_y;
				
				pts_x[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
				pts_y[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));
			}
			
			// create a spline
			tk::spline s;
			
			// set (x, y) points to the spline
			s.set_points(pts_x, pts_y);
			
			// Define the actual (x, y) points we will use for the planner
			vector<double> next_x_vals;
			vector<double> next_y_vals;
			
			// Start with all of the previous path points from last time
			for(int i=0; i<previous_path_x.size(); ++i)
			{
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}
			
			// Calculate how to break up spline points so that we travel at our desired reference velocity
			double tgt_x = (double)spacing;
			double tgt_y = s(tgt_x);
			double tgt_dist = sqrt((tgt_x)*(tgt_x) + (tgt_y)*(tgt_y));
			
			double x_add_on = 0;
			
			for(int i=1; i<=60-previous_path_x.size(); ++i)
			{
				double N = (tgt_dist / (.02*curr_ref_vel/2.24));
				double x_point = x_add_on + (tgt_x)/N;
				double y_point = s(x_point);
				
				x_add_on = x_point;
				
				double x_ref = x_point;
				double y_ref = y_point;
				
				// rotate back to normal after rotating it earlier
				x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
				y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));
				
				x_point += ref_x;
				y_point += ref_y;
				
				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}
			
			// END
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
