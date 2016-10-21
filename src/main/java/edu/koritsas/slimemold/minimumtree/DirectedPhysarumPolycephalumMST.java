package edu.koritsas.slimemold.minimumtree;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.*;
import org.geotools.graph.structure.basic.BasicDirectedNode;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by ilias on 21/10/2016.
 */
public abstract class DirectedPhysarumPolycephalumMST extends PhysarumPolycephalumMST{
    public DirectedPhysarumPolycephalumMST(DirectedGraph graph, Node sourceNode, List<Node> sinkNodes, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNode, sinkNodes, Io, γ, numberOfIterations);

    }
    @Override
    protected RealVector createConstantsVector() {
        List<Node> nodes = getAllNodes(graph);

        double[] constants =new double[nodes.size()];
        for (int i=0;i<nodes.size();i++){
            Node n=nodes.get(i);
            if (n.equals(sourceNode)) {
                constants[i] = Io;
            }else{
                constants[i]=-Io/(sinkNodesList.size());
            }

        }


        return  new ArrayRealVector(constants);
    }


    @Override
    protected double calculateCoefficient(Node n1, Node n2) {
        BasicDirectedNode node1 = (BasicDirectedNode) n1;
        BasicDirectedNode node2 = (BasicDirectedNode) n2;
        List<Edge> ndf =node1.getInEdges();

        double c=0;
        Edge edge =node1.getEdge(node2);
        if(node1.getOutEdges().contains(edge)){

            c=super.calculateCoefficient(node1,node2);
        }else if(node1.getInEdges().contains(edge)){
            c=-FastMath.abs(super.calculateCoefficient(node1,node2));
        }else{
            c=0;
        }


        return c;

    }

    @Override
    protected double calculateTubeFlux(Edge e) {
        DirectedEdge edge = (DirectedEdge) e;
        double p1 =pressureMap.get(edge.getInNode());
        double p2 =pressureMap.get(edge.getOutNode());

        double D =conductivityMap.get(edge);
        double w = getEdgeCost(e);

        double Q=(D/w)*(p1-p2);
        if(Q<0){
            Q=0;
        }


        return Q;

    }
    @Override
    public double calculateTubeConductivity(Edge e) {
        DirectedEdge edge= (DirectedEdge) e;
        double D = conductivityMap.get(edge);


        double Q = fluxMap.get(e);


        double L = getEdgeCost(e);

       // double p1 = pressureMap.get(edge.getInNode());
       // double p2 = pressureMap.get(edge.getOutNode());

       // double ps = pressureMap.get(sourceNode);
      //  double pe = pressureMap.get(sinkNode);


        double fQ = FastMath.pow(FastMath.abs(Q),γ);


        double newD = fQ - 0.1 * D;





        return newD;
    }
}
