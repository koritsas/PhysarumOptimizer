package edu.koritsas.slimemold.mstree;

import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumLangrarianCSP;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

import java.util.Iterator;
import java.util.List;

/**
 * Created by ilias on 7/11/2016.
 */
public abstract class PhysarumPolycephalumLagrarianCSPT extends PhysarumPolycephalumLangrarianCSP {
    public PhysarumPolycephalumLagrarianCSPT(Graph graph, Node sourceNode, List<Node> sinkNodes, double absoluteThreshold, double relativeThreshold, int numberOfIterations, double λmax, double step) {
        super(graph, sourceNode, null, absoluteThreshold, relativeThreshold, numberOfIterations, λmax, step);
        this.sinkNodesList=sinkNodes;

    }



    @Override
    protected RealVector createConstantsVector() {
        List<Node> allNodes =getAllNodes(graph);
        double[] constants = new double[allNodes.size()];

        for (int i = 0; i < constants.length; i++) {

            Node n1 = allNodes.get(i);
            if (n1.equals(sourceNode)) {

                constants[i] = Io;
            } else if (sinkNodesList.contains(n1)) {
                constants[i] = -Io/sinkNodesList.size();
            } else {

                constants[i] = 0;

            }
        }


        return new ArrayRealVector(constants);


    }


    @Override
    public double calculateTubeConductivity(Edge e) {
        double D = conductivityMap.get(e);

        double Q = currentFluxMap.get(e);

        double fQ= FastMath.abs(Q);
        double newD = fQ - 0.45 * D;


        return newD;
    }
}
