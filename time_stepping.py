import numpy as np
import matplotlib.pyplot as plt
import math

def f(y):
	return -2 * y

def analytical(t):
	return 0.5 * math.exp(-2*t)

def explicitEulerEq(y, dt):
	return y[-1] + dt * f(y[-1])

def explicitEuler(y0, dt, tn):
	y = np.array([y0])
	t = 0.0
	while t + dt < tn:
		y = np.append(y, explicitEulerEq(y, dt))
		t += dt
	y = np.append(y, explicitEulerEq(y, (tn - t)))
	return y[-1]

def leapfrogEq(y, dt):
	return y[-2] + 2.0 * dt * f(y[-1])

def leapfrog(y0, dt, tn):
	y = np.array([y0])
	y = np.append(y, explicitEulerEq(y, dt)) #leapfrog need tn-1 and tn, do not have information yet, then uses exp euler for that
	t = 0.0
	while t + dt < tn:
		y = np.append(y, leapfrogEq(y, dt))
		t += dt
	y = np.append(y, leapfrogEq(y, (tn - t)))
	return y[-1]

def rungeKutta4Eq(y, dt):
	k1 = f(y[-1])
	k2 = f(y[-1] + dt / 2.0 * k1)
	k3 = f(y[-1] + dt / 2.0 * k2)
	k4 = f(y[-1] + dt * k3)
	return y[-1] + dt / 6.0 * (k1 + 2*k2 + 2*k3 + k4)

def rungeKutta4(y0, dt, tn):
	y = np.array([y0])
	t = 0.0
	while t + dt < tn:
		y = np.append(y, rungeKutta4Eq(y, dt))
		t += dt
	y = np.append(y, rungeKutta4Eq(y, (tn - t)))
	return y[-1]

def error(n1, n2):
	return math.fabs(n1 - n2)

def getDtList():
	step = 0.01
	return np.arange(step, 0.64, step)
	# return np.array([0.64,0.32,0.16,0.08,0.04,0.02,0.01,0.005,0.0025,0.00125,0.000625,0.0003125,0.00015625, 7.8125e-5,3.90625e-5,1.953125e-5,9.765625e-6,4.8828125e-6])

def main():
	dt_list = getDtList()
	y0 = 0.5
	tn = 2.0
	analytical_answer = analytical(tn)
	exp_euler_error = np.array([])
	leapfrog_error = np.array([])
	rk4_error = np.array([])
	for dt in dt_list:
		exp_euler_error = np.append(exp_euler_error, error(explicitEuler(y0, dt, tn), analytical_answer))
		leapfrog_error = np.append(leapfrog_error, error(leapfrog(y0, dt, tn), analytical_answer))
		rk4_error = np.append(rk4_error, error(rungeKutta4(y0, dt, tn), analytical_answer))

	plt.subplot(2, 2, 1)
	plt.plot(dt_list, exp_euler_error)
	plt.title("exp euler")

	plt.subplot(2, 2, 2)
	plt.plot(dt_list, leapfrog_error)
	plt.title("leapfrog")

	plt.subplot(2, 2, 3)
	plt.plot(dt_list, rk4_error)
	plt.title("rk4")

	# plt.legend(['exp euler', 'leapfrog', 'rk4'])
	plt.show()

main()