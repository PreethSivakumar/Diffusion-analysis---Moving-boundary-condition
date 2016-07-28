
/*	Given the concentration profiles and interface position at timestep j and estimates for the values at j+1,
	updates the estimate for interface position at timestep j+1 */
int new_interface_planar(two_phase *tp, double time_step, double pass_nmb);
int new_interface_spherical(two_phase *tp, double time_step, double pass_nmb);

/*	Given the concentration profiles and interface position at timestep j and estimates for the interface at j+1,
	finds concentrations at j+1 in phase to the left of the boundary	*/
void new_concentration_left_planar(two_phase *tp, double time_step);
void new_concentration_left_spherical(two_phase *tp, double time_step);

/*	Given the concentration profiles and interface position at timestep j and estimates for the interface at j+1,
	finds concentrations at j+1 in phase to the right of the boundary	*/
void new_concentration_right_planar(two_phase *tp, double time_step);
void new_concentration_right_spherical(two_phase *tp, double time_step);		


/*	Takes system from time step j to timestep j+1 */
int	take_step_planar(two_phase *tp, double dt, double tolerance);
int	take_step_spherical(two_phase *tp, double dt, double tolerance);
