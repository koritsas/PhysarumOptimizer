package edu.koritsas.slimemold.shapefile;

import org.geotools.graph.structure.Graph;

/**
 * Created by koritsas on 15/12/2016.
 */
public abstract class Constraint {
    private Graph mGraph;
    private final double mMinValue;
    private final double mMaxValue;


    public Constraint(Graph solutionGraph, double constraintMinValue, double constraintMaxValue){
        this.mGraph=solutionGraph;
        this.mMinValue=constraintMinValue;
        this.mMaxValue=constraintMaxValue;
    }
    public double getmMinValue() {
        return mMinValue;
    }

    public double getmMaxValue() {
        return mMaxValue;
    }

    abstract double getConstraintValue();

    public boolean isViolated(){
        boolean violated=false;
        if (getConstraintValue()<=mMinValue || getConstraintValue() >=mMaxValue){
            violated=true;
        }
        return violated;
    }
}
