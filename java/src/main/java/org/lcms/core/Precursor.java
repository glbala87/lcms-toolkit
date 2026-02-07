package org.lcms.core;

/**
 * Precursor ion information for MS/MS spectra.
 */
public class Precursor {
    private double mz;
    private double intensity;
    private int charge;
    private double isolationWindowLower;
    private double isolationWindowUpper;
    private ActivationMethod activationMethod;
    private double collisionEnergy;

    public Precursor() {
        this.activationMethod = ActivationMethod.UNKNOWN;
    }

    public Precursor(double mz, int charge) {
        this();
        this.mz = mz;
        this.charge = charge;
    }

    // Getters and setters
    public double getMz() { return mz; }
    public void setMz(double mz) { this.mz = mz; }

    public double getIntensity() { return intensity; }
    public void setIntensity(double intensity) { this.intensity = intensity; }

    public int getCharge() { return charge; }
    public void setCharge(int charge) { this.charge = charge; }

    public boolean hasCharge() { return charge != 0; }

    public double getIsolationWindowLower() { return isolationWindowLower; }
    public void setIsolationWindowLower(double lower) { this.isolationWindowLower = lower; }

    public double getIsolationWindowUpper() { return isolationWindowUpper; }
    public void setIsolationWindowUpper(double upper) { this.isolationWindowUpper = upper; }

    public double getIsolationWindowWidth() {
        return isolationWindowUpper - isolationWindowLower;
    }

    public ActivationMethod getActivationMethod() { return activationMethod; }
    public void setActivationMethod(ActivationMethod method) { this.activationMethod = method; }

    public double getCollisionEnergy() { return collisionEnergy; }
    public void setCollisionEnergy(double energy) { this.collisionEnergy = energy; }

    /**
     * Calculate neutral mass from m/z and charge.
     */
    public double getNeutralMass() {
        if (charge == 0) return 0.0;
        double protonMass = 1.007276;
        return (mz - protonMass) * Math.abs(charge);
    }

    @Override
    public String toString() {
        return String.format("Precursor(mz=%.4f, charge=%d)", mz, charge);
    }
}
