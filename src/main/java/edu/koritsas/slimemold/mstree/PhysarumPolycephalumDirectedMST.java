package edu.koritsas.slimemold.mstree;

import edu.koritsas.slimemold.AbstractDirectedPhysarumPolycephalum;
import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.geotools.graph.structure.*;

import java.util.List;

/**
 * Created by ilias on 21/10/2016.
 */
public abstract class PhysarumPolycephalumDirectedMST extends AbstractDirectedPhysarumPolycephalum {
    public PhysarumPolycephalumDirectedMST(DirectedGraph graph, Node sourceNode, List<Node> sinkNodes, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNode, null, Io, γ, numberOfIterations);
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
        DirectedEdge edge= (DirectedEdge) e;
        double D = conductivityMap.get(edge);


        double Q = fluxMap.get(e);


        double L = getEdgeCost(e);

        //double fQ = Math.pow(Math.abs(Q), γ) / (1 + Math.pow(Math.abs(Q), γ));
         double fQ=Math.abs(Q);

         double newD = fQ - 0.2 * D;



        return newD;
    }
}
