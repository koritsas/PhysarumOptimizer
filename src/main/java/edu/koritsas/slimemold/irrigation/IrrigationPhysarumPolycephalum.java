package edu.koritsas.slimemold.irrigation;

import edu.koritsas.slimemold.mstree.PhysarumPolycephalumLagrarianCSPT;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.graph.traverse.standard.DijkstraIterator;
import org.opengis.feature.simple.SimpleFeature;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by koritsas on 9/1/2017.
 */
public abstract class IrrigationPhysarumPolycephalum extends PhysarumPolycephalumLagrarianCSPT {
    protected HashMap<Constraint,Double> constraintHashMap = new HashMap<>();
    public IrrigationPhysarumPolycephalum(Graph graph, Node sourceNode, List<Node> sinkNodes, double absoluteThreshold, double relativeThreshold, int numberOfIterations, double λmax, double step) {
        super(graph, sourceNode, sinkNodes, absoluteThreshold, relativeThreshold, numberOfIterations, λmax, step);
    }
    public void createConstraints(){

        DijkstraIterator.EdgeWeighter weigter = new DijkstraIterator.EdgeWeighter() {
            @Override
            public double getWeight(Edge e) {
                return 0;
            }
        };

        DijkstraShortestPathFinder pf = new DijkstraShortestPathFinder(graph,sourceNode,weigter);
        pf.calculate();

        for (Node n:sinkNodesList){
         constraintHashMap.putIfAbsent(new Constraint(pf.getPath(n)),λ);

        }

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
            selfC=selfC+ FastMath.abs(-D/w);
        }

        return selfC;
    }


    @Override
    public double calculateTubeConductivity(Edge e) {
        double D = conductivityMap.get(e);

        double Q = currentFluxMap.get(e);

        double fQ= FastMath.abs(Q);
        double newD = fQ - 0.8 * D;


        return newD;
    }

    @Override
    public void execute() {
        double startingTime=System.currentTimeMillis();
        createConstraints();
        initializeMaps(graph);
        logger.info("Starting iterations...");
        while (!converged) {
            //for (int i = 0; i <200 ; i++) {

            List<Edge> edges = new ArrayList<>(graph.getEdges());
            edges.stream().forEach(edge -> L.put(edge, L.get(edge) + λ * getEdgeConstraintValue(edge)));
            for (Edge e:edges){
                Iterator<Constraint> it=constraintHashMap.keySet().iterator();
                while (it.hasNext()){
                    Constraint c=it.next();
                    if(c.containsEdge(e)){
                        L.put(e,L.get(e)+constraintHashMap.get(c)*getEdgeConstraintValue(e));
                    }
                }

            }


            while (!converged) {




                logger.info("Current iteration: " + iteration);
                List<Node> allNodes = getAllNodes(graph);


                List<Node> allButSinkNodes = getAllButSink(graph);


                RealVector constants = createConstantsVector();

                RealMatrix coefficients = createCoefficientsMatrix();


                // DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

                DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();

                RealVector solution = solver.solve(constants);

                redefinePressures(allButSinkNodes, solution);

                redefineFlows(graph);

                redefineDiameters(graph);


                converged = checkConvergenceCriterion(previousFluxMap, currentFluxMap, iteration, m_absoluteThreshold, m_relativeThreshold);


                previousFluxMap = new HashMap<>(currentFluxMap);
                iteration++;
                if (iteration==numberOfIterations){
                    break;
                }

            }


            if (iteration==numberOfIterations){
                break;
            }
            Graph newGraph=getGraph();



           // converged=false;
            List<Edge> edges1 =new ArrayList<>(newGraph.getEdges());

            //boolean conv =pathViolatesConstraints(newGraph);

            Iterator<Constraint> iter = constraintHashMap.keySet().iterator();
            while (iter.hasNext()){
                Constraint con =iter.next();
                if (con.isViolated(getGraph())){
                    converged=false;
                    constraintHashMap.put(con,constraintHashMap.get(con)+step);
                }

            }


         initializeMaps(graph);
        }
        double endingTime=System.currentTimeMillis();

        logger.info("Execution Time: "+ ((endingTime-startingTime))/1000);
    }
}
