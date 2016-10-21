package edu.koritsas.slimemold.mstree;

import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by ilias on 21/10/2016.
 */
public abstract class PhysarumPolycephalumMST extends AbstractPhysarumPolycephalum {

    public PhysarumPolycephalumMST(Graph graph, Node sourceNode, List<Node> sinkNodes, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNode,null, Io, γ, numberOfIterations);
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

        double Q = fluxMap.get(e);
        double L=getEdgeCost(e);



        //double fQ = Math.pow(Math.abs(Q), γ) / (1 + Math.pow(Math.abs(Q), γ));
        double fQ=Math.abs(Q);

        double newD = fQ - 0.4 * D;



        return newD;
    }

}
