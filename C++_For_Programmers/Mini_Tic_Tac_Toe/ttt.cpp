// Main file for Tic Tac Toe game.
#include "ttt.hpp"
#include <iostream>


int main() {

    Board tic_tac_toe;

    while (!tic_tac_toe.complete && tic_tac_toe.rounds < 9){

        // Decide who's turn it is.
        char player;
        if (tic_tac_toe.rounds % 2 == 0) {
            player = 'X';
        }
        else {
            player = 'O';
        }

        // Prints the board.
        tic_tac_toe.display_board(tic_tac_toe.board);

        // Collect the users move and checks that the move is legal. 
        int move = tic_tac_toe.collect_input(player);

        // Updates the board.
        tic_tac_toe.update_board(tic_tac_toe.board, move, player);

        // Checks for a winner.
        if (tic_tac_toe.check_for_winner(tic_tac_toe.board)) {
            std::cout << "PLAYER " << player << " WINS!!\n";
            tic_tac_toe.display_board(tic_tac_toe.board);
            tic_tac_toe.complete = true;
            break;

        }

        // Checks to see if we have a draw.
        if (tic_tac_toe.rounds == 8) {
            std::cout << "Game is a draw!\n";
            tic_tac_toe.display_board(tic_tac_toe.board);
            tic_tac_toe.complete = true;
            break;
        }

        tic_tac_toe.rounds++;
    }

    return 0;
}
