# Outage Probability Analysis of OTFS-Based Downlink LEO Satellite Communication with UAV Cooperative Relaying

**B.Tech Senior Design Project (SDP)**
**Jacalyn Rena Karri | VIT-AP University, Amaravati | 2026**
**Guide: Dr. M. Mahesh | School of Electronics Engineering**

---

## Overview

Low Earth Orbit (LEO) satellite systems experience severe Doppler shifts due to high orbital speeds, which degrade conventional OFDM-based communication systems. This project analytically investigates **outage probability** for an **OTFS-modulated** downlink LEO satellite communication system enhanced by **UAV cooperative relaying**.

The key idea: OTFS modulation operates in the **delay-Doppler domain**, converting rapidly varying satellite channels into nearly time-invariant representations. A UAV relay operating in **Decode-and-Forward (DF) mode** provides an additional transmission path, significantly improving communication reliability.

---

## System Architecture

```
         LEO Satellite (Source Node)
         Ns = 10 transmit antennas
              |              |
   Direct Link|              | Satellite → UAV Link
  (Shadowed-  |              | (Shadowed-Rician Fading)
   Rician)    |              |
              |           UAV Relay
              |        (Decode-and-Forward)
              |              |
              |              | UAV → Destination Link
              |              | (Nakagami-m Fading)
              |______________|
                     |
            Destination Node
            (Ground User, Single Antenna)
            MRC combining of both paths
```

**Two transmission paths:**
- Direct: LEO Satellite → Destination Node
- Relay: LEO Satellite → UAV Relay → Destination Node

---

## OTFS Signal Processing Chain

```
Transmitter:                    Channel:              Receiver:
DD Domain  → ISFFT → Heisenberg → SR + Nakagami-m → Wigner → SFFT → FD-LE
x[k,l]       (TF    (Time         Fading +                           Equalizer
              Domain) Domain)     Doppler                            → x̂[k,l]
```

- **ISFFT** — Inverse Symplectic Finite Fourier Transform
- **Heisenberg Transform** — converts TF domain to time domain
- **SFFT** — Symplectic FFT at receiver
- **FD-LE** — Frequency Domain Linear Equalization (suppresses Doppler interference)

---

## Channel Models

| Link | Fading Model | Reason |
|------|-------------|--------|
| LEO Satellite → Destination | Shadowed-Rician | Captures both LoS and shadowing in satellite channels |
| LEO Satellite → UAV | Shadowed-Rician | Same satellite link characteristics |
| UAV → Destination | Nakagami-m | Models terrestrial relay fading conditions |

---

## Outage Probability Framework

Outage occurs when achievable mutual information falls below target rate Rb:

```
P_out = P(I_OTFS < Rb)
```

The outage variable is defined as:
```
X_OTFS = L_sd * Z_sd + ξ * L_rd * Z_rd
P_out = P(X_OTFS < η),   η = σ²(2^(2Rb) - 1) / Ps
```

### Statistical Approximation Steps

| Step | Operation |
|------|-----------|
| 1 | Shadowed-Rician fading → compute mean and variance using hypergeometric functions |
| 2 | Moment-matching → approximate χ_sd ~ Gamma(α, β) |
| 3 | Reciprocal variable → Z_sd ~ Inverse-Gamma(α_z, β_z) |
| 4 | Inverse-Gamma property → 1/Z_sd ~ Gamma |
| 5 | Final closed-form → P_out = γ(α_OTFS, η·β_OTFS) / Γ(α_OTFS) |

**Final closed-form outage probability:**
```
P_out = γ(α_OTFS, η·β_OTFS) / Γ(α_OTFS)
```
where γ(·,·) is the lower incomplete Gamma function implemented via MATLAB's `gammainc()`.

---

## Simulation Parameters

| Parameter | Description | Value |
|-----------|-------------|-------|
| N | Delay bins | 4 |
| M | Doppler bins | 4 |
| Ns | Satellite transmit antennas | 10 |
| b_sd | Shadowed-Rician parameter | 0.251 |
| m_sd | Shadowed-Rician parameter | 0.251 |
| Ω_sd | Average channel power (SD link) | 0.279 |
| Ω_rd | Average channel power (RD link) | 1 |
| σ² | Noise variance | 1 |
| Rb | Target transmission rate | 1 bps/Hz |
| m_rd | Nakagami relay fading parameter | 2, 3, 4 |
| ξ | Relay cooperation factor | 0.3, 0.5 |
| SNR | Transmit SNR range | 0 – 25 dB |

---

## Key Results

- Outage probability **decreases with increasing transmit SNR** — confirming that higher transmit power improves reliability
- **UAV cooperative relaying significantly outperforms** direct transmission without relay
- Larger Nakagami-m parameter (m_rd: 2 → 3 → 4) → **faster outage reduction** due to less severe relay fading
- Larger relay cooperation factor ξ → **stronger relay contribution** → better outage performance
- OTFS provides **robust performance under severe Doppler conditions** compared to conventional OFDM

---

## File Structure

```
OTFS-LEO-UAV-Outage/
├── analytical_outage_otfs.m      # Main analytical outage probability computation
├── relay_outage_simulation.m     # Relay fading case simulation
├── otfs_expression_test.m        # Expression validation and testing
└── README.md                     # This file
```

---

## How to Run

1. Open MATLAB
2. Run `analytical_outage_otfs.m`
3. Output: Outage probability vs SNR curves for different relay fading parameters (m_rd = 2, 3, 4) and cooperation factors (ξ = 0.3, 0.5)

**Requirements:** MATLAB with Statistics and Signal Processing Toolboxes

---

## Concepts Covered

- OTFS modulation and delay-Doppler domain signal processing
- LEO satellite channel modeling (Shadowed-Rician fading)
- UAV cooperative relaying (Decode-and-Forward mode)
- Frequency Domain Linear Equalization (FD-LE)
- Statistical moment-matching approximation (Gamma and Inverse-Gamma distributions)
- Closed-form outage probability derivation using incomplete Gamma functions
- MATLAB analytical framework implementation

---

## Reference

Jia Shi, Junfan Hu, Yang Yue, Xuan Xue, Wei Liang, and Zan Li,
**"Outage Probability for OTFS Based Downlink LEO Satellite Communication,"**
*IEEE Transactions on Vehicular Technology*, vol. 71, no. 3, pp. 3355–3360, March 2022.

---

## Author

**Jacalyn Rena Karri**
B.Tech Electronics and Communication Engineering (Minor: Automotive Design)
VIT-AP University, Amaravati — 2026
GitHub: [github.com/Jacalynrena14](https://github.com/Jacalynrena14)
