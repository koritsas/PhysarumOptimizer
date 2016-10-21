package edu.shortestpath.test;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;

import java.util.*;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Created by ilias on 5/10/2016.
 */
public abstract class AbstractPhysarumPolycephalum {
    protected Graph graph;
    protected Node sourceNode;
    protected Node sinkNode;
    protected List<Node> sinkNodesList;
    protected double Io;
    protected double γ;
    protected  int numberOfIterations;
    protected HashMap<Edge,Double> fluxMap ;
    protected HashMap<Edge,Double> conductivityMap;
    protected HashMap<Node,Double> pressureMap;

    protected static final Logger logger = Logger.getLogger(AbstractPhysarumPolycephalum.class.getName());

    protected int iteration=0;

    protected HashMap<Edge, XYSeries> flowSeriesMap = new HashMap<Edge, XYSeries>();
    protected HashMap<Edge, XYSeries> conductivitySeriesMap =new HashMap<Edge, XYSeries>();
    protected HashMap<Node, XYSeries> pressureSeriesMap =new HashMap<Node, XYSeries>();

    public HashMap<Edge, XYSeries> getFlowSeriesMap() {
        return flowSeriesMap;
    }

    public void setFlowSeriesMap(HashMap<Edge, XYSeries> flowSeriesMap) {
        this.flowSeriesMap = flowSeriesMap;
    }

    public HashMap<Edge, XYSeries> getConductivitySeriesMap() {
        return conductivitySeriesMap;
    }

    public void setConductivitySeriesMap(HashMap<Edge, XYSeries> conductivitySeriesMap) {
        this.conductivitySeriesMap = conductivitySeriesMap;
    }

    /**
     *
     * @return the graph
     */
    public Graph getGraph() {
        return graph;
    }



    /**
     * Executes the algorithm
     */
    public void execute() {

        initializeMaps(graph);

        logger.info("Starting iterations...");

        for (int i = 0; i <numberOfIterations ; i++) {

            logger.info("Current iteration: "+i);
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


            iteration++;
        }

        eliminateEdges(graph);
    }


    /**
     *
     * @return the coefficients matrix A for the system Ax=B;
     */


    protected RealMatrix createCoefficientsMatrix() {

        List<Node> allNodes=getAllNodes(graph);
        List<Node> allButSink =getAllButSink(graph);

        double[][] matrix = new double[allNodes.size()][allButSink.size()];

        for (int i = 0; i < allNodes.size(); i++) {
            Node n1 = allNodes.get(i);
            for (int j = 0; j < allButSink.size(); j++) {
                Node n2 = allButSink.get(j);

                if(n1.equals(n2)){
                    matrix[i][j]=calculateSelfCoefficient(n1);
                }else{
                        matrix[i][j]=calculateCoefficient(n1,n2);
                }


            }
        }


        return new Array2DRowRealMatrix(matrix);
    }

    /**
     *
     * @return the constants vector B for the system Ax=B
     */

    protected RealVector createConstantsVector() {
        List<Node> allNodes =getAllNodes(graph);
        double[] constants = new double[allNodes.size()];
        for (int i = 0; i < constants.length; i++) {

            Node n1 = allNodes.get(i);
            if (n1.equals(sourceNode)) {

                    constants[i] = Io - calculateCoefficient(n1,sinkNode) * pressureMap.get(sinkNode);

            } else if (n1.equals(sinkNode)) {
                constants[i] = -Io - calculateSelfCoefficient(n1) *  pressureMap.get(sinkNode);
            } else {

                    constants[i] = -calculateCoefficient(n1,sinkNode) *  pressureMap.get(sinkNode);


            }
        }


        return new ArrayRealVector(constants);
    }
    /**
     * Initializes the XYSeries of the diagrams
     * @param graph
     */
    protected void initializeMaps(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        fluxMap=new HashMap<Edge, Double>(graph.getEdges().size());
        conductivityMap=new HashMap<Edge,Double>(graph.getEdges().size());
        for (Edge e:edges){
            fluxMap.putIfAbsent(e,0.0);
            conductivityMap.putIfAbsent(e,1.0);
        }

        Collection<Node> nodes =graph.getNodes();
        pressureMap=new HashMap<Node,Double>(nodes.size());
        for (Node n:nodes){
            pressureMap.putIfAbsent(n,0.0);

        }




    }



    /**
     * @param solution Redefines pressures on edges from the solution of the system
     */
    protected void redefinePressures(List<Node> allButSink, RealVector solution) {
        for (int i = 0; i < allButSink.size(); i++) {

            Node n = allButSink.get(i);
            //noinspection Since15
            pressureSeriesMap.putIfAbsent(n,new XYSeries("("+n.getID()+")"));
            pressureMap.put(n,solution.getEntry(i));
            pressureSeriesMap.get(n).add(iteration,solution.getEntry(i));
        }


    }

    /**
     *
     * @param graph
     * redefines the flows for the edges of the graph
     */
    @SuppressWarnings("Since15")
    protected void redefineFlows(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        for(Edge e:edges){
            flowSeriesMap.putIfAbsent(e,new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double Q =calculateTubeFlux(e);
            fluxMap.put(e,Q);
            flowSeriesMap.get(e).add(iteration,Q);
        }

    }

    /**
     *
     * @return the cost of the solution
     */
    public double getSolutionCost() {

        Collection<Edge> edges = graph.getEdges();
        double cost = edges.stream().collect(Collectors.summingDouble(value -> getEdgeCost(value)));
        return  cost;
    }

    /**
     *
     * @param graph eliminates the edges of the graph that have conductivity near zero
     */
    protected void eliminateEdges(Graph graph){

        graph.getEdges().removeIf(o -> FastMath.abs(fluxMap.get(o))<10E-10);


    }


    /**
     *
     * @param graph
     * redefines the flows for the edges of the graph
     */
    protected void redefineDiameters(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        for(Edge e:edges){
            //noinspection Since15
            conductivitySeriesMap.putIfAbsent(e, new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double D =calculateTubeConductivity(e);
            conductivityMap.put(e,D);
            conductivitySeriesMap.get(e).add(iteration,D);
        }

    }

    /**
     *
     * @param e
     * @return the flux of the tube
     */
    protected double calculateTubeFlux(Edge e){

        double p1 =pressureMap.get(e.getNodeA());
        double p2 =pressureMap.get(e.getNodeB());

        double D =conductivityMap.get(e);
        double w = getEdgeCost(e);

        double Q=(D/w)*(p1-p2);

        return Q;
    }

    /**
     *
     * @param e
     * @return the new conductivity of the edge
     */
    public abstract double calculateTubeConductivity(Edge e);



    /**
     * show the pressure diagram
     */
    public void showPressureMap(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );

/*
        Collection<Node> nodes =graph.getNodes();
        for(Node n:nodes){
            dataset.addSeries(pressureSeriesMap.get(n));
        }
    */
        Iterator it= pressureSeriesMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }
        JFreeChart flowChart = ChartFactory.createXYLineChart("Time-Pressure", "Time", "H", dataset, PlotOrientation.VERTICAL,true,true,true);
        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");

        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);


    }
    public void showConductivityMap(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );
        Iterator it= conductivitySeriesMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }

        JFreeChart flowChart =ChartFactory.createXYLineChart("Time-Conductivity", "Time", "D", dataset, PlotOrientation.VERTICAL,true,false,false);
        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);


    }


    public void showFlowDiagram(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );

        Iterator it= flowSeriesMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }
        JFreeChart flowChart =ChartFactory.createXYLineChart("Time-Flow", "Time", "Q", dataset,PlotOrientation.VERTICAL,true,true,true);

        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    protected List<Node> getAllButSink(Graph graph) {




        return getAllNodes(graph).stream().filter(node -> node!=sinkNode).collect(Collectors.toList());
    }


    /**
     * @param graph
     * @return Returns a list with all the nodes in the graph
     */

    protected List<Node> getAllNodes(Graph graph) {


        List<Node> nodeList = new ArrayList<>(graph.getNodes());


        return nodeList;

    }

    /**
     * @param n
     * @return Returns the self-coefficient
     */

    public double calculateSelfCoefficient(Node n) {
        double selfC = 0;
        List<Edge> edges = n.getEdges();

        Iterator<Node> it =n.getRelated();

        while(it.hasNext()){
            selfC = selfC + FastMath.abs(calculateCoefficient(n,it.next()));
        }


        return selfC;
    }

    protected double calculateCoefficient(Node n1, Node n2){
        Edge e=n1.getEdge(n2);
        double coeff=0;
        if (e!=null){
            double D =conductivityMap.get(e);
            double w = getEdgeCost(e);
            coeff=-D/w;
        }

        return coeff;
    }

    /**
     *
     * @param e
     * @return the weight of the edge e
     */
    public abstract double getEdgeCost(Edge e);







}
