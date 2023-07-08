// Source file for defining classes and functions for Tic Tac Toe game.
#include "ttt.hpp"
#include <iostream>
#include <string>

// Board Class Implementations

// Constructor.
Board::Board() 
    : board{' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '}, complete(false), rounds(0) {}

// Collect the user move.
int Board::collect_input(char player) {
    
    int move;
    std::cout << "Player " << player << " Enter a value from 1-9:  ";
    std::cin >> move;
    while(!this->check_legal_move(this->board, move)){
        std::cout << "Invalid input, please select another:  ";
        std::cin >> move;
    }
    std::cout << "\n";
    return move;
}

// Check that the move is legal (used in collect_input()).
bool Board::check_legal_move(char *board, int move) {
    // Decrement the move (input as 1-9 but array is indexed 0-8)
    int move_index = move - 1;
    // Checks that the move suggested by the user is legal.
    if ((board[move_index] == ' ') && (move > 0) && (move < 10)) {
        return true;
    }
    return false;
}

// Apply the input to the board / update the board with the user input.
void Board::update_board(char *board, int move, char player) {
    // Updates the board with the user input.
    int move_index = move - 1;
    board[move_index] = player;
}

// Display the board on the console.
void Board::display_board(char *board) {
    std::cout << "\n";
    std::cout << "     |     |    \n";
    std::cout << "  " << board[0] << "  |  " << board[1] << "  |  " << board[2] << "\n";
    std::cout << "____ | ___ | ____\n";
    std::cout << "     |     |    \n";
    std::cout << "  " << board[3] << "  |  " << board[4] << "  |  " << board[5] << "\n";
    std::cout << "____ | ___ | ____\n";
    std::cout << "     |     |    \n";
    std::cout << "  " << board[6] << "  |  " << board[7] << "  |  " << board[8] << "\n";
    std::cout << "     |     |    \n";
    std::cout << "\n";
}

// Check for a winner.
bool Board::check_for_winner(char *board) {

    // Check horizontal.
    if ((board[0] != ' ') && (board[0] == board[1]) && (board[1] == board[2])) {
        return true;
    }
    else if ((board[3] != ' ') && (board[3] == board[4]) && (board[4] == board[5])) {
        return true;
    }
    else if ((board[6] != ' ') && (board[6] == board[7]) && (board[7] == board[8])) {
        return true;
    }

    // Check vertical.
    if ((board[0] != ' ') && (board[0] == board[3]) && (board[3] == board[6])) {
        return true;
    }
    else if ((board[1] != ' ') && (board[1] == board[4]) && (board[4] == board[7])) {
        return true;
    }
    else if ((board[2] != ' ') && (board[2] == board[5]) && (board[5] == board[8])) {
        return true;
    }

    // Check diagonal. 
    if ((board[0] != ' ') && (board[0] == board[4]) && (board[4] == board[8])) {
        return true;
    }
    else if ((board[2] != ' ') && (board[2] == board[4]) && (board[4] == board[6])) {
        return true;
    }

    return false;
}

