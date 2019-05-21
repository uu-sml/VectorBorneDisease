#!/usr/bin/python
"""Processing of the simulation data"""
import json
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats


# Parameters
burnin = 500
input_folder = "output/"
file_pgibbs = input_folder + "yap_dengue_pgibbs_2048.json"
file_csmc = input_folder + "yap_dengue_csmc_2048.json"

output_folder = "data/"
traces_file = output_folder + "vbd_traces.csv"
density_file = output_folder + "vbd_density.csv"
autocorr_file = output_folder + "vbd_autocorr.csv"
histogram_file = output_folder + "vbd_histogram.csv"

# Load pgibbs data
samples_pgibbs = []
rho_csmc = []
with open(file_pgibbs, "r") as f:
    data = json.load(f)
    for d in data[burnin:-1]:
        samples_pgibbs.append(d["θ"]["ρ"])

# Load csmc data
w_csmc = []
with open(file_csmc, "r") as f:
    data = json.load(f)
    for d in data[burnin:-1]:
        rho_csmc.append(d["θ"]["ρ"])  # Careful!
        w_csmc.append(d["lweight"])

# Normalize weights
w_csmc = np.exp(w_csmc - np.max(w_csmc))
w_csmc /= np.sum(w_csmc)

# Find posterior as mixture of beta
rho = np.linspace(0.1, 0.35, 100)
posterior_csmc = np.zeros(rho.size)
samples_csmc = []
for (r, w) in zip(rho_csmc, w_csmc):
    alpha = r["α"]
    beta = r["β"]
    samples_csmc.append(scipy.stats.beta.rvs(alpha, beta))
    posterior_csmc += w*scipy.stats.beta.pdf(rho, alpha, beta)


def corr(x, n_lags):
    N = len(x)
    autocorr = np.zeros(n_lags)

    mu = np.mean(x)

    for k in range(n_lags):
        autocorr[k] = np.dot(x[k:] - mu, x[:N-k] - mu)/(N-k)
    autocorr /= autocorr[0]
    return autocorr


n_lags = 100
corr_csmc = corr(samples_csmc, n_lags)
corr_pgibbs = corr(samples_pgibbs, n_lags)


# PLOTTING
def plot():

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    ax1.plot(rho, posterior_csmc)
    ax1.hist(samples_pgibbs, bins=15, density=True)
    ax1.set_ylabel("Density")
    ax1.set_xlabel("Reporting rate")
    ax1.legend(["mPG", "PG"])

    ax2.plot(samples_csmc)
    ax2.plot(samples_pgibbs)
    ax2.set_xlabel("Sample")
    ax2.set_ylabel("Reporting rate")
    ax2.legend(["mPG", "PG"])

    ax3.plot(corr(samples_csmc, n_lags))
    ax3.plot(corr(samples_pgibbs, n_lags))
    ax3.set_xlabel("Lag")
    ax3.set_ylabel("Reporting rate")
    ax3.legend(["mPG", "PG"])

    plt.show()


def save():

    bins, locations = np.histogram(samples_pgibbs, bins=27, density=True)

    with open(traces_file, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(["sample", "PG", "mPG"])
        writer.writerows(zip(range(len(samples_pgibbs)),
                             samples_pgibbs, samples_csmc))

    with open(density_file, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(["rho", "density"])
        writer.writerows(zip(rho, posterior_csmc))

    with open(autocorr_file, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(["lag", "PG", "mPG"])
        writer.writerows(zip(range(n_lags), corr_pgibbs, corr_csmc))

    with open(histogram_file, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(["rho", "density"])
        writer.writerows(zip(locations[0:-1], bins))


if __name__ == "__main__":
    plot()
    save()
