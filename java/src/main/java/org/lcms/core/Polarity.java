package org.lcms.core;

/**
 * Ion polarity mode.
 */
public enum Polarity {
    UNKNOWN(0),
    POSITIVE(1),
    NEGATIVE(-1);

    private final int value;

    Polarity(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public static Polarity fromValue(int value) {
        for (Polarity p : values()) {
            if (p.value == value) return p;
        }
        return UNKNOWN;
    }
}
