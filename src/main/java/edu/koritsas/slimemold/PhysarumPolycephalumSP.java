package edu.koritsas.slimemold;

import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

/**
 * Created by ilias on 5/10/2016.
 */

public abstract class PhysarumPolycephalumSP extends AbstractPhysarumPolycephalum{
    public PhysarumPolycephalumSP(Graph graph, Node sourceNode,Node sinkNode,int numberOfIterations){
        this.graph=graph;
        this.sourceNode=sourceNode;
        this.sinkNode=sinkNode;
        this.numberOfIterations=numberOfIterations;

    }

    @Override
    public double calculateTubeConductivity(Edge e) {
        return 0;
    }
}
