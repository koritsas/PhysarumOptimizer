package edu.koritsas.slimemold;

import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumSP;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

/**
 * Created by ilias on 9/11/2016.
 */
public abstract class TestSlime extends PhysarumPolycephalumSP{
    public TestSlime(Graph graph, Node sourceNode, Node sinkNode, double absoluteThreshold, double relativeThreshold, int numberOfIterations) {
        super(graph, sourceNode, sinkNode, absoluteThreshold, relativeThreshold, numberOfIterations);
    }

    @Override
    protected double calculateCoefficient(Node n1, Node n2) {
        Edge e=n1.getEdge(n2);
        double coeff=0;
        if (e!=null){
            double D =conductivityMap.get(e);
            double w = getEdgeCost(e);

            coeff=-D/w;
        }

        return coeff;
    }

    public abstract double getEdgeConstraint(Edge edge);
}
