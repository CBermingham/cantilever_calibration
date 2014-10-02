def Sho_power_response(omega, omega_R, Q, A, B):
	return A + (B * omega_R^4) / ((omega^2 - omega_R^2)^2 + ((omega^2 * omega_R^2) / Q^2))