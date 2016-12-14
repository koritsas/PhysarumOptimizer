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
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

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
    protected double calculateCoefficient(Node n1, Node n2){

        //Edge e=n1.getEdge(n2);


       /* if (e!=null){
            double D =conductivityMap.get(e);
            double w = getEdgeCost(e);
            coeff=-D/w;
        }
*/
        List<Edge> edges = n1.getEdges(n2);
        double coeff=0;

        if (edges!=null){
            Iterator<Edge> edgeIterator =edges.iterator();
            while (edgeIterator.hasNext()){
                Edge edge=edgeIterator.next();
                double D =conductivityMap.get(edge);
                double w = L.get(edge);
                coeff=coeff-D/w;
            }

        }
        return coeff;
    }

    @Override
    public double calculateSelfCoefficient(Node n) {
        double selfC = 0;
        List<Edge> edges = n.getEdges();


       /* Iterator<Node> it =n.getRelated();

        while(it.hasNext()){
            selfC = selfC + FastMath.abs(calculateCoefficient(n,it.next()));
        }
*/
        for (Edge edge:edges){
            double D =conductivityMap.get(edge);
            double w = L.get(edge) ;
            selfC=selfC+FastMath.abs(-D/w);
        }

        return selfC;
    }


    @Override
    public double calculateTubeConductivity(Edge e) {
        double D = conductivityMap.get(e);

        double Q = currentFluxMap.get(e);

        double fQ= FastMath.abs(Q);
        double newD = fQ - 0.4 * D;


        return newD;
    }
}
