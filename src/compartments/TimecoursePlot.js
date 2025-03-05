import React, { useState, useEffect } from "react";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";
import PropTypes from "prop-types";

const TimecoursePlot = ({ gene, dataUrl }) => {
    const [chartData, setChartData] = useState([]);

    useEffect(() => {
        const fetchData = async () => {
            try {
                const response = await fetch(dataUrl);
                const jsonData = await response.json();
                if (jsonData.genes && jsonData.genes[gene]) {
                    const rawData = jsonData.genes[gene];
                    const regions = ["cortex", "interface", "medulla"];
                    const aggregatedData = {};

                    regions.forEach(region => {
                        if (rawData[region]) {
                            rawData[region].forEach(([time, value]) => {
                                if (!aggregatedData[time]) {
                                    aggregatedData[time] = { time };
                                }
                                if (!aggregatedData[time][region]) {
                                    aggregatedData[time][region] = [];
                                }
                                aggregatedData[time][region].push(value);
                            });
                        }
                    });

                    // Compute averages
                    const formattedData = Object.keys(aggregatedData).map(time => {
                        const result = { time: Number(time) };
                        regions.forEach(region => {
                            if (aggregatedData[time][region]) {
                                result[region] = aggregatedData[time][region].reduce((a, b) => a + b, 0) / aggregatedData[time][region].length;
                            }
                        });
                        return result;
                    }).sort((a, b) => a.time - b.time); // Ensure sorting by time

                    setChartData(formattedData);
                }
            } catch (error) {
                console.error("Error fetching timecourse data:", error);
            }
        };

        fetchData();
    }, [gene, dataUrl]);

    return (
        <ResponsiveContainer width="100%" height={400}>
            <LineChart data={chartData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="time" label={{ value: "Time (hours)", position: "insideBottom", offset: -5 }} domain={[0, 48]} ticks={[0, 12, 24, 48]} />
                <YAxis label={{ value: "Average Expression", angle: -90, position: "insideLeft" }} />
                <Tooltip />
                <Legend />
                <Line type="monotone" dataKey="cortex" stroke="#8884d8" strokeWidth={2} />
                <Line type="monotone" dataKey="interface" stroke="#82ca9d" strokeWidth={2} />
                <Line type="monotone" dataKey="medulla" stroke="#ff7300" strokeWidth={2} />
            </LineChart>
        </ResponsiveContainer>
    );
};

TimecoursePlot.propTypes = {
    gene: PropTypes.string.isRequired,
    dataUrl: PropTypes.string.isRequired,
};

export default TimecoursePlot;
